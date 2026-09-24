version 1.0

import "../../structs/Structs.wdl"

task FilterVcfForHmmIBD {

    meta {
        description: "Pre-filter a VCF for IBD inference: mask low-depth genotypes (FMT/DP < min_depth), optionally preserve the original allele-frequency annotations, recompute AN/AC/AF (optionally per population), and drop alleles/sites that are no longer represented. Emits a compressed BCF ready for hmmibd-rs --from-bcf."

        tool:          "bcftools"
        tool_version:  "1.22"
        tool_url:      "https://www.htslib.org/"
        tool_citation: "Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008."

        author: "Jonn Smith"

        outputs: {
            filtered_bcf:       "Filtered, AN/AC/AF-recomputed callset as a compressed BCF",
            filtered_bcf_index: "CSI index for filtered_bcf"
        }
    }

    parameter_meta {
        input_vcf:             "VCF/BCF to filter prior to IBD inference. (required)"
        prefix:                "Basename for the output BCF. (required)"
        min_depth:             "Genotypes with FORMAT/DP below this value are set to missing (bcftools filter -S .). (default: 5)"
        keep_original_af:      "If true, rename the existing INFO annotations via rename_annots_tsv before recomputing AN/AC/AF, so the original frequencies are preserved under new tags. (default: false)"
        variant_types:         "Which variant types to keep (bcftools view -v): 'snps', 'indels', or 'both'. SNPs are the standard hmmIBD marker set; Pf indels are error-prone. (default: snps)"
        biallelic_only:        "Restrict to biallelic sites (bcftools view -m2 -M2). Recommended true: multiallelic sites (especially indels) can exceed hmmibd-rs --max-all and crash it until that is patched. (default: true)"
        split_multiallelics:   "Split multiallelic records into biallelic ones (bcftools norm -m-any) so SNP alleles at multiallelic/spanning-deletion sites are recovered rather than dropped. Runs after the type pre-select to stay fast. (default: true)"
        max_variants:          "Optional cap on the number of variants kept. When >0 and the filtered callset exceeds it, sites are thinned evenly across the genome down to at most this many (not truncated to the first N). (default: 0 = no limit)"
        mask_gq0_genotypes:    "Also set GQ0 genotypes to missing (in addition to the FORMAT/DP < min_depth mask). Older GATK (pre-4.6.0.0 GenotypeGVCFs; also GnarlyGenotyper) emitted no-/low-confidence hom-refs as 0/0 with GQ=0 instead of ./. (GATK issue #7792, fixed in PR #8741). Turning this on reproduces that fix downstream — the DP mask alone misses GQ0 calls whose DP >= min_depth. NOTE: this is the conservative choice — it drops ALL GQ0 hom-refs, including any that are genuinely well-covered ref; use a GVCF cross-reference if you need to keep those. (default: false)"
        populations_file:      "Optional sample-to-population file passed to `bcftools +fill-tags -S`; when given, AN/AC/AF are computed per population. When omitted, tags are computed across all samples. (default: none)"
        rename_annots_tsv:     "Required when keep_original_af=true: two-column TSV of old-name<TAB>new-name passed to `bcftools annotate --rename-annots`. (default: none)"
        extra_args:            "Additional `bcftools view` filters appended to the FINAL view (after AN/AC/AF are recomputed, so INFO-based expressions like MAF work). Whitespace-separated tokens; write expressions without internal spaces, e.g. \"-e MAF<0.01\" or \"-i F_MISSING<0.1\". Note the final view already applies -i 'MAX(INFO/AC) > 0'; supply extra exclusions with -e (a second -i would override the built-in one). (default: empty)"
        runtime_attr_override: "Override the default runtime attributes. (default: none)"
    }

    input {
        File input_vcf
        String prefix

        Int min_depth = 5
        Boolean keep_original_af = false
        String variant_types = "snps"
        Boolean biallelic_only = true
        Boolean split_multiallelics = true
        Int max_variants = 0
        Boolean mask_gq0_genotypes = false

        File? populations_file
        File? rename_annots_tsv

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(input_vcf, "GB"))

    # Genotype mask: always drop low-depth genotypes; optionally also drop GQ0 genotypes
    # (see mask_gq0_genotypes / GATK issue #7792). bcftools sets matching genotypes to ./. (-S .).
    String gt_mask_expr = if mask_gq0_genotypes then "FMT/DP < " + min_depth + " | FMT/GQ = 0" else "FMT/DP < " + min_depth
    # FORMAT fields to KEEP in the up-front strip. GQ is normally dropped, but must survive when
    # mask_gq0_genotypes needs it in gt_mask_expr below (otherwise bcftools errors: GQ not in header).
    String fmt_keep = if mask_gq0_genotypes then "^FORMAT/GT,FORMAT/AD,FORMAT/DP,FORMAT/GQ" else "^FORMAT/GT,FORMAT/AD,FORMAT/DP"

    command <<<
        set -euxo pipefail

        # ---- Resource detection (required preamble) ----
        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        # Reserve 1 GB for OS + container overhead.
        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        # Per-thread RAM (used by samtools sort -m, etc.)
        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        # Java heap (used by GATK/Picard). Same as USABLE_RAM_GB; alias for clarity.
        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"
        # ---- end preamble ----

        # Variant-type and biallelic selection, kept as separate flags so the type can be
        # pre-selected BEFORE `norm` (below). hmmibd-rs can technically use indels (it reads
        # allele indices, not allele strings), but Pf indels are error-prone, so SNPs are the
        # default. Keeping multiallelic sites (biallelic_only=false) can exceed hmmibd-rs
        # --max-all and crash it until the fork patch lands.
        # Assign to a variable first so shellcheck sees a variable (not a constant) in `case`,
        # and keep the flag sets as arrays so their expansion is quote-safe (shellcheck-clean).
        VTYPE="~{variant_types}"
        case "${VTYPE}" in
            snps)   TYPE_FLAGS=(-v snps) ;;
            indels) TYPE_FLAGS=(-v indels) ;;
            both)   TYPE_FLAGS=() ;;
            *) echo "ERROR: variant_types must be one of: snps, indels, both" >&2 ; exit 1 ;;
        esac
        BIALLELIC_FLAGS=()
        if [[ "~{biallelic_only}" == "true" ]] ; then BIALLELIC_FLAGS=(-m2 -M2) ; fi

        # User-supplied extra `bcftools view` filters (e.g. "-e MAF<0.01" or "-i F_MISSING<0.1").
        # Word-split into an array so multiple tokens pass through cleanly under set -u and quoted
        # expansion. Tokens are whitespace-separated, so write expressions without internal spaces
        # (bcftools accepts e.g. -e MAF<0.01). Applied in the FINAL view, after AN/AC/AF are
        # recomputed, so INFO-based expressions work. Empty by default.
        read -r -a EXTRA_ARGS <<< "~{extra_args}"

        # Split multiallelic records into biallelic ones, so the alleles at a multiallelic or
        # spanning-deletion site are recovered instead of the whole site being discarded by the
        # biallelic filter. This runs for EVERY variant_types setting (snps, indels, or both) --
        # it is not gated to SNPs. It runs after the type pre-select purely for speed: selecting
        # a type first means norm only splits records of the type(s) you keep (so snps mode never
        # touches the raw ~250-allele indel sites). With variant_types=both nothing is pre-dropped,
        # so norm splits everything -- correct, just slower on multiallelic-heavy cohorts. Runs
        # after the up-front annotation strip, so norm never sees PL or unneeded INFO. `cat` is a
        # no-op passthrough when disabled.
        NORM_STEP=(cat)
        if [[ "~{split_multiallelics}" == "true" ]] ; then NORM_STEP=(bcftools norm -m-any -Ou) ; fi

        # keeping the original allele frequencies requires a rename map
        if [[ "~{keep_original_af}" == "true" && -z "~{rename_annots_tsv}" ]] ; then
            echo "ERROR: keep_original_af=true requires rename_annots_tsv to be provided." >&2
            exit 1
        fi

        # Mask low-depth genotypes, (optionally) preserve original AF via rename, recompute
        # AN/AC/AF (optionally per population), then trim now-unrepresented alt alleles/sites.
        if [[ "~{keep_original_af}" == "true" ]] ; then
            bcftools annotate -x '~{fmt_keep}' -Ou ~{input_vcf} \
              | bcftools filter -S . -e "~{gt_mask_expr}" -Ou \
              | bcftools view "${TYPE_FLAGS[@]}" -Ou \
              | "${NORM_STEP[@]}" \
              | bcftools annotate --rename-annots ~{rename_annots_tsv} -Ou \
              | bcftools +fill-tags -Ou -- ~{"-S " + populations_file} -t AN,AC,AF \
              | bcftools view "${BIALLELIC_FLAGS[@]}" "${TYPE_FLAGS[@]}" --trim-alt-alleles -i 'MAX(INFO/AC) > 0' -Ob --threads "${NUM_CPUS}" -o ~{prefix}.filtered.bcf
        else
            # First step: strip to only the fields anything downstream uses. hmmibd-rs reads
            # FORMAT/GT (first-ploidy), FORMAT/AD (dominant-allele), and the FILTER column;
            # our DP-mask needs FORMAT/DP. Everything else — all INFO, and FORMAT PL (Number=G,
            # quadratic in alleles) / GQ / phasing — is dropped. This is much smaller and ~3x
            # faster on big cohorts, and it removes GATK per-allele INFO (e.g. HAPCOMP) with
            # value counts that disagree with the ALT count and would abort --trim-alt-alleles.
            bcftools annotate -x 'INFO,~{fmt_keep}' -Ou ~{input_vcf} \
              | bcftools filter -S . -e "~{gt_mask_expr}" -Ou \
              | bcftools view "${TYPE_FLAGS[@]}" -Ou \
              | "${NORM_STEP[@]}" \
              | bcftools +fill-tags -Ou -- ~{"-S " + populations_file} -t AN,AC,AF \
              | bcftools view "${BIALLELIC_FLAGS[@]}" "${TYPE_FLAGS[@]}" --trim-alt-alleles -i 'MAX(INFO/AC) > 0' -Ob --threads "${NUM_CPUS}" -o ~{prefix}.filtered.bcf
        fi

        # Apply user-supplied extra filters as a SEPARATE view: `bcftools view` accepts only one
        # -i/-e expression and the view above already uses -i 'MAX(INFO/AC) > 0', so a user -i/-e
        # (or any other view-level filter) must run as its own pass here, on the recomputed callset.
        if [[ ${#EXTRA_ARGS[@]} -gt 0 ]] ; then
            bcftools view "${EXTRA_ARGS[@]}" -Ob --threads "${NUM_CPUS}" ~{prefix}.filtered.bcf -o ~{prefix}.filtered.extra.bcf
            mv ~{prefix}.filtered.extra.bcf ~{prefix}.filtered.bcf
        fi

        bcftools index --threads "${NUM_CPUS}" ~{prefix}.filtered.bcf

        # Optional cap on the number of variants (off when max_variants <= 0). When the
        # filtered callset has more sites than max_variants, thin it evenly across the
        # genome (keep every ceil(total / max_variants)-th site) so the retained variants
        # stay genome-wide, instead of truncating to the first N.
        if [[ ~{max_variants} -gt 0 ]] ; then
            TOTAL=$(bcftools index -n ~{prefix}.filtered.bcf)
            echo "variant cap: max_variants=~{max_variants}, filtered site count=${TOTAL}"
            if [[ "${TOTAL}" -gt ~{max_variants} ]] ; then
                STRIDE=$(( (TOTAL + ~{max_variants} - 1) / ~{max_variants} ))
                bcftools view ~{prefix}.filtered.bcf \
                  | awk -v s="${STRIDE}" -v m=~{max_variants} '
                        /^#/ { print; next }
                        { n++; if (((n - 1) % s) == 0 && k < m) { print; k++ } }' \
                  | bcftools view -Ob --threads "${NUM_CPUS}" -o ~{prefix}.filtered.capped.bcf
                mv ~{prefix}.filtered.capped.bcf ~{prefix}.filtered.bcf
                bcftools index -f --threads "${NUM_CPUS}" ~{prefix}.filtered.bcf
                echo "variant cap: retained $(bcftools index -n ~{prefix}.filtered.bcf) variants (stride ${STRIDE})"
            fi
        fi
    >>>

    output {
        File filtered_bcf       = "~{prefix}.filtered.bcf"
        File filtered_bcf_index = "~{prefix}.filtered.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task HmmIBDrs {

    meta {
        description: "Run hmmibd-rs (the Rust reimplementation of hmmIBD) directly on a BCF to infer identity-by-descent segments and per-pair IBD fractions. Every hmmibd-rs flag is exposed as an input with its upstream default, including the output-mode switches (suppress_frac, bcf_to_bin_file, bcf_to_bin_file_by_chromosome) that change which files are produced. Outputs are glob-backed so each is present only when the chosen flags actually write it: a normal run yields ibd_segments + ibd_fraction; suppress_frac drops ibd_fraction; a bcf-to-bin mode skips the HMM and yields binary_genotypes instead."

        tool:          "hmmibd-rs"
        tool_version:  "0.1.5"
        tool_url:      "https://github.com/bguo068/hmmibd-rs"
        tool_citation: "Guo B, Schaffner SF, Taylor AR, et al. hmmibd-rs. Malar J. 2026;25:110."

        author: "Jonn Smith"

        outputs: {
            ibd_segments:     "Per sample-pair IBD/non-IBD segments with assigned state and variant counts (<prefix>.hmm.txt); empty when a bcf-to-bin mode skips the HMM",
            ibd_fraction:     "Per sample-pair summary; key column fract_sites_IBD (<prefix>.hmm_fract.txt); empty when suppress_frac=true or a bcf-to-bin mode is used",
            binary_genotypes: "Binary genotype file(s) (<prefix>*.bin); non-empty only under bcf_to_bin_file or bcf_to_bin_file_by_chromosome"
        }
    }

    parameter_meta {
        input_bcf:             "Genotype input for hmmibd-rs. With from_bcf=true (default): a BCF/VCF read via --from-bcf. With from_bcf=false: an hmmIBD text genotype table (-i), e.g. from BcfToSampleTable. (required)"
        prefix:                "Output prefix; produces <prefix>.hmm.txt and <prefix>.hmm_fract.txt. (required)"
        from_bcf:              "Read the input as BCF/VCF (--from-bcf, applies --bcf-read-mode). Set false to read a pre-built hmmIBD text genotype table instead. (default: true)"

        data_file2:            "Optional second-population genotypes (-I/--data-file2). (default: none)"
        freq_file1:            "Optional allele-frequency file for population 1 (-f/--freq-file1); computed from data when omitted. (default: none)"
        freq_file2:            "Optional allele-frequency file for population 2 (-F/--freq-file2). (default: none)"
        bad_samples_file:      "Optional list of sample IDs to exclude (-b/--bad-file). (default: none)"
        good_pairs_file:       "Optional list of sample pairs to analyze (-g/--good-file). (default: none)"
        bcf_filter_config:     "Optional TOML BCF filter config (--bcf-filter-config); a default is generated when omitted. (default: none)"
        genome:                "Optional genome/recombination-map spec (--genome); mutually exclusive with rec_rate, and takes precedence when set. (default: none)"

        bcf_read_mode:         "How to read genotypes from the BCF (--bcf-read-mode): dominant-allele | first-ploidy | each-ploidy. (default: dominant-allele)"

        max_iter:              "Max EM iterations (-m/--max-iter). (default: 5)"
        k_rec_max:             "Cap on inferred number of generations (-n/--k-rec-max). (default: none; hmmibd-rs uses Inf = no cap)"
        eps:                   "Genotyping error rate (--eps). (default: 0.001)"
        min_inform:            "Minimum number of informative sites for a pair (--min-inform). (default: 10)"
        min_discord:           "Minimum discordance fraction to consider a pair (--min-discord). (default: 0.0)"
        max_discord:           "Maximum discordance fraction to consider a pair (--max-discord). (default: 1.0)"
        min_snp_sep:           "Minimum bp separation between SNPs (--min-snp-sep). (default: 5)"
        fit_thresh_dpi:        "Convergence threshold on pi (--fit-thresh-dpi). (default: 0.001)"
        fit_thresh_dk:         "Convergence threshold on k (--fit-thresh-dk). (default: 0.01)"
        fit_thresh_drelk:      "Convergence threshold on relative k (--fit-thresh-drelk). (default: 0.001)"
        rec_rate:              "Constant recombination rate per generation per bp (-r/--rec-rate); ignored when genome is set. (default: 0.00000074)"

        max_all:               "Maximum number of unique alleles per site (--max-all). (default: 8)"
        buffer_size_segments:  "Write buffer size in bytes for the segments output (--buffer-size-segments). (default: none; hmmibd-rs uses 8Kb)"
        buffer_size_frac:      "Write buffer size in bytes for the fraction output (--buffer-size-frac). (default: none; hmmibd-rs uses 8Kb)"

        filt_min_seg_cm:       "Drop output segments shorter than this many cM (--filt-min-seg-cm). (default: none)"
        filt_max_tmrca:        "Drop pairs with k_rec above this threshold (--filt-max-tmrca). (default: none)"
        filt_ibd_only:         "Omit non-IBD segments from the output (--filt-ibd-only). (default: false)"

        num_threads:           "Number of threads; 0 = all CPUs (--num-threads). (default: 0)"
        par_mode:              "Parallelization mode: 0 = small sample sets, 1 = large sample sets (--par-mode). (default: 0)"
        par_chunk_size:        "Sample-pairs per chunk (par-mode 0) or samples per chunk (par-mode 1) (--par-chunk-size). (default: 120)"

        suppress_frac:                 "Suppress the per-pair fraction output (--suppress-frac); leaves ibd_fraction empty. (default: false)"
        bcf_to_bin_file:               "Convert the BCF to a single binary genotype file and skip the HMM (--bcf-to-bin-file); yields binary_genotypes, leaves ibd_segments/ibd_fraction empty. (default: false)"
        bcf_to_bin_file_by_chromosome: "Convert the BCF to one binary genotype file per chromosome and skip the HMM (--bcf-to-bin-file-by-chromosome); yields binary_genotypes, leaves ibd_segments/ibd_fraction empty. (default: false)"

        extra_args:            "Additional command-line args appended verbatim to the hmmibd-rs invocation. (default: empty)"
        runtime_attr_override: "Override the default runtime attributes. (default: none)"
    }

    input {
        File input_bcf
        String prefix

        File? data_file2
        File? freq_file1
        File? freq_file2
        File? bad_samples_file
        File? good_pairs_file
        File? bcf_filter_config
        File? genome

        Boolean from_bcf = true
        String bcf_read_mode = "dominant-allele"

        Int max_iter = 5
        Float? k_rec_max
        Float eps = 0.001
        Int min_inform = 10
        Float min_discord = 0.0
        Float max_discord = 1.0
        Int min_snp_sep = 5
        Float fit_thresh_dpi = 0.001
        Float fit_thresh_dk = 0.01
        Float fit_thresh_drelk = 0.001
        String rec_rate = "0.00000074"

        Int max_all = 8
        Int? buffer_size_segments
        Int? buffer_size_frac

        Float? filt_min_seg_cm
        Float? filt_max_tmrca
        Boolean filt_ibd_only = false

        Int num_threads = 0
        Int par_mode = 0
        Int par_chunk_size = 120

        Boolean suppress_frac = false
        Boolean bcf_to_bin_file = false
        Boolean bcf_to_bin_file_by_chromosome = false

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(input_bcf, "GB"))

    # --genome and -r/--rec-rate are mutually exclusive; genome wins when supplied.
    String recombination_arg = if defined(genome) then "--genome " + select_first([genome]) else "-r " + rec_rate
    # --bcf-read-mode only applies when reading a BCF/VCF; omit it in text-table mode.
    String read_mode_arg = if from_bcf then "--bcf-read-mode " + bcf_read_mode else ""

    command <<<
        set -euxo pipefail

        # ---- Resource detection (required preamble) ----
        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        # Reserve 1 GB for OS + container overhead.
        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        # Per-thread RAM (used by samtools sort -m, etc.)
        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        # Java heap (used by GATK/Picard). Same as USABLE_RAM_GB; alias for clarity.
        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"
        # ---- end preamble ----

        hmmibd-rs \
            ~{true="--from-bcf" false="" from_bcf} \
            -i ~{input_bcf} \
            -o ~{prefix} \
            ~{read_mode_arg} \
            ~{"-I " + data_file2} \
            ~{"-f " + freq_file1} \
            ~{"-F " + freq_file2} \
            ~{"-b " + bad_samples_file} \
            ~{"-g " + good_pairs_file} \
            ~{"--bcf-filter-config " + bcf_filter_config} \
            -m ~{max_iter} \
            ~{"-n " + k_rec_max} \
            --eps ~{eps} \
            --min-inform ~{min_inform} \
            --min-discord ~{min_discord} \
            --max-discord ~{max_discord} \
            --min-snp-sep ~{min_snp_sep} \
            --fit-thresh-dpi ~{fit_thresh_dpi} \
            --fit-thresh-dk ~{fit_thresh_dk} \
            --fit-thresh-drelk ~{fit_thresh_drelk} \
            ~{recombination_arg} \
            --max-all ~{max_all} \
            ~{"--buffer-size-segments " + buffer_size_segments} \
            ~{"--buffer-size-frac " + buffer_size_frac} \
            ~{"--filt-min-seg-cm " + filt_min_seg_cm} \
            ~{"--filt-max-tmrca " + filt_max_tmrca} \
            ~{true="--filt-ibd-only" false="" filt_ibd_only} \
            --num-threads ~{num_threads} \
            --par-mode ~{par_mode} \
            --par-chunk-size ~{par_chunk_size} \
            ~{true="--suppress-frac" false="" suppress_frac} \
            ~{true="--bcf-to-bin-file" false="" bcf_to_bin_file} \
            ~{true="--bcf-to-bin-file-by-chromosome" false="" bcf_to_bin_file_by_chromosome} \
            ~{extra_args}
    >>>

    output {
        # Glob so each output exists only when the chosen flags actually write it.
        Array[File] ibd_segments     = glob("~{prefix}.hmm.txt")
        Array[File] ibd_fraction     = glob("~{prefix}.hmm_fract.txt")
        Array[File] binary_genotypes = glob("~{prefix}*.bin")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-hmmibd-rs:0.1.5"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task BcfToVcf {

    meta {
        description: "Re-emit a BCF as a bgzipped, tabix-indexed VCF (bcftools view -Oz). Convenience conversion so the filtration workflow can hand back a VCF."

        tool:          "bcftools"
        tool_version:  "1.22"
        tool_url:      "https://www.htslib.org/"

        author: "Jonn Smith"

        outputs: {
            filtered_vcf:       "The input BCF re-emitted as a bgzipped VCF",
            filtered_vcf_index: "Tabix (.tbi) index for filtered_vcf"
        }
    }

    parameter_meta {
        input_bcf:             "BCF (or VCF) to re-emit as bgzipped VCF. (required)"
        prefix:                "Basename for the output VCF. (required)"
        runtime_attr_override: "Override the default runtime attributes. (default: none)"
    }

    input {
        File input_bcf
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(input_bcf, "GB"))

    command <<<
        set -euxo pipefail
        NUM_CPUS=$(grep -c '^processor' /proc/cpuinfo)

        bcftools view "~{input_bcf}" -Oz --threads "${NUM_CPUS}" -o "~{prefix}.filtered.vcf.gz"
        bcftools index -t --threads "${NUM_CPUS}" "~{prefix}.filtered.vcf.gz"
    >>>

    output {
        File filtered_vcf       = "~{prefix}.filtered.vcf.gz"
        File filtered_vcf_index = "~{prefix}.filtered.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task BcfToSampleTable {

    meta {
        description: "Convert a bi-allelic BCF/VCF to the hmmIBD text genotype table: a tab-delimited matrix with columns chrom(int) pos then one call per sample, alleles coded 0/1/... and -1 for missing, using the first allele of each GT (matches hmmibd-rs --bcf-read-mode first-ploidy). Also emits a matching bi-allelic allele-frequency file (same sites and order) computed from the matrix, so it can be supplied to both hmmIBD and hmmibd-rs. Non-numeric contigs (MIT/API) are dropped because hmmIBD requires integer chromosomes."

        tool:          "bcftools"
        tool_version:  "1.22"
        tool_url:      "https://www.htslib.org/"

        author: "Jonn Smith"

        outputs: {
            sample_gt_table:   "hmmIBD text genotype table (<prefix>.sample_gt.txt); hmmibd-rs input in text mode (from_bcf=false)",
            sample_freq_table: "Bi-allelic allele-frequency file (<prefix>.sample_freq.txt), same sites/order as the genotype table"
        }
    }

    parameter_meta {
        input_bcf:             "Bi-allelic, SNP-filtered BCF/VCF (as produced by FilterVcfForHmmIBD) to convert. (required)"
        prefix:                "Basename for the output tables. (required)"
        runtime_attr_override: "Override the default runtime attributes. (default: none)"
    }

    input {
        File input_bcf
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(input_bcf, "GB"))

    command <<<
        set -euxo pipefail

        GT="~{prefix}.sample_gt.txt"
        FRQ="~{prefix}.sample_freq.txt"

        # header: chrom  pos  <sample1> <sample2> ...
        { printf 'chrom\tpos'; bcftools query -l "~{input_bcf}" | awk '{printf "\t%s", $0}'; printf '\n'; } > "${GT}"

        # body: integer chrom, pos, first-allele call per sample (-1 missing)
        bcftools query -f '%CHROM\t%POS[\t%GT]\n' "~{input_bcf}" | awk 'BEGIN{FS=OFS="\t"}
        {
            c=$1; sub(/^Pf3D7_/,"",c); sub(/_v[0-9]+$/,"",c)   # Pf3D7_01_v3 -> 01
            if (c !~ /^[0-9]+$/) next                          # drop non-numeric contigs
            printf "%d\t%d", c+0, $2
            for (i=3;i<=NF;i++){ n=split($i,a,/[\/|]/); g=a[1]; if (g=="."||g=="") g=-1; printf "\t%s", g }
            printf "\n"
        }' >> "${GT}"

        # bi-allelic allele frequencies from the matrix (identical sites and order)
        awk 'NR==1{next}
        {
            c0=0; c1=0; n=0
            for (i=3;i<=NF;i++){ if($i=="0"){c0++;n++} else if($i=="1"){c1++;n++} }
            if (n==0){ f0=1; f1=0 } else { f0=c0/n; f1=c1/n }
            printf "%s\t%s\t%.6f\t%.6f\n", $1, $2, f0, f1
        }' "${GT}" > "${FRQ}"

        echo "variants: $(( $(wc -l < "${GT}") - 1 )); samples: $(( $(head -1 "${GT}" | awk '{print NF}') - 2 ))"
    >>>

    output {
        File sample_gt_table   = "~{prefix}.sample_gt.txt"
        File sample_freq_table = "~{prefix}.sample_freq.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task StitchHmmIBDTables {

    meta {
        description: "Combine many per-region hmmIBD genotype tables (one per input VCF, as produced by BcfToSampleTable) into a SINGLE hmmIBD table + matching allele-frequency file for hmmibd-rs. This is the scalable combine step: each huge input VCF is first reduced to its compact integer genotype table, and only these small tables are stitched here — a merged multi-hundred-GB BCF is never materialized. Assumes every input table carries the SAME samples in the SAME column order (i.e. the input VCFs are disjoint regions of ONE joint call, split across files); the sample headers are verified identical and the task aborts if they differ. Rows (sites) from all tables are concatenated and coordinate-sorted by integer chrom then pos (hmmibd-rs assumes sorted input); the allele-frequency file is recomputed on the COMBINED matrix (valid because every table holds the full sample set) so no cross-file reconciliation is needed. An optional max_variants cap is applied to the combined, genome-wide callset (thinned by even stride, not truncated). Input regions are expected to be non-overlapping; same-position rows (e.g. a SNP and an indel from variant_types=both) are left in place and handled downstream by hmmibd-rs --min-snp-sep."

        tool:          "coreutils"
        tool_url:      "https://www.gnu.org/software/coreutils/"

        author: "Jonn Smith"

        outputs: {
            sample_gt_table:   "Combined hmmIBD text genotype table (<prefix>.sample_gt.txt): all input sites, coordinate-sorted, one first-allele call per sample (-1 missing)",
            sample_freq_table: "Combined bi-allelic allele-frequency file (<prefix>.sample_freq.txt), same sites/order as the genotype table, recomputed across all samples"
        }
    }

    parameter_meta {
        gt_tables:             "hmmIBD genotype tables to combine, one per input VCF (from BcfToSampleTable / FilterVcfForHmmIBD output_format='hmmibd_table'). All must share the same samples in the same column order. (required)"
        prefix:                "Basename for the combined output tables. (required)"
        max_variants:          "Optional cap on the number of variants in the COMBINED callset. When >0 and the combined site count exceeds it, sites are thinned evenly across the genome down to at most this many (not truncated to the first N). Applied genome-wide here (not per input file), so the cap is meaningful across the merged callset. (default: 0 = no limit)"
        runtime_attr_override: "Override the default runtime attributes. (default: none)"
    }

    input {
        Array[File] gt_tables
        String prefix
        Int max_variants = 0
        RuntimeAttr? runtime_attr_override
    }

    # Sort is external (disk-backed); size for input + sort temp + output.
    Int disk_size = 20 + ceil(10.0 * size(gt_tables, "GB"))

    command <<<
        set -euxo pipefail

        GT="~{prefix}.sample_gt.txt"
        FRQ="~{prefix}.sample_freq.txt"

        TABLES=(~{sep=' ' gt_tables})

        # Header (chrom pos <samples...>) comes from the first table; every other table MUST have
        # the identical header. The tables are disjoint regions of ONE joint call, so the sample
        # set and column order are the same across files; a mismatch means the inputs are not the
        # same samples and stitching would misalign genotype columns -> abort loudly.
        head -n1 "${TABLES[0]}" > header.txt
        for t in "${TABLES[@]}"; do
            if ! head -n1 "$t" | cmp -s - header.txt ; then
                echo "ERROR: sample header mismatch in $t; all input tables must share the same samples in the same column order (disjoint regions of one joint call)." >&2
                exit 1
            fi
        done

        # Concatenate all bodies (drop each header) and coordinate-sort by integer chrom then pos.
        # External sort (-T .) keeps temp files on the task disk, so this scales past RAM.
        for t in "${TABLES[@]}"; do tail -n +2 "$t"; done \
            | sort -T . -k1,1n -k2,2n > body_sorted.txt

        # Write header, then the (optionally thinned) sorted body.
        cat header.txt > "${GT}"
        if [[ ~{max_variants} -gt 0 ]]; then
            TOTAL=$(wc -l < body_sorted.txt)
            echo "combined variant cap: max_variants=~{max_variants}, combined site count=${TOTAL}"
            if [[ "${TOTAL}" -gt ~{max_variants} ]]; then
                # Even genome-wide stride over the sorted callset (keep every STRIDE-th site),
                # matching FilterVcfForHmmIBD's per-file cap logic but applied to the merged set.
                STRIDE=$(( (TOTAL + ~{max_variants} - 1) / ~{max_variants} ))
                awk -v s="${STRIDE}" -v m=~{max_variants} '{ if (((NR-1) % s) == 0 && k < m) { print; k++ } }' body_sorted.txt >> "${GT}"
                echo "combined variant cap: retained $(( $(wc -l < "${GT}") - 1 )) variants (stride ${STRIDE})"
            else
                cat body_sorted.txt >> "${GT}"
            fi
        else
            cat body_sorted.txt >> "${GT}"
        fi
        rm -f body_sorted.txt

        # Recompute bi-allelic allele frequencies on the COMBINED matrix (explicit freq file,
        # same sites/order as GT). Valid because every input table carries the full sample set,
        # so a per-site frequency over the merged rows equals concatenating per-file frequencies.
        awk 'NR==1{next}
        {
            c0=0; c1=0; n=0
            for (i=3;i<=NF;i++){ if($i=="0"){c0++;n++} else if($i=="1"){c1++;n++} }
            if (n==0){ f0=1; f1=0 } else { f0=c0/n; f1=c1/n }
            printf "%s\t%s\t%.6f\t%.6f\n", $1, $2, f0, f1
        }' "${GT}" > "${FRQ}"

        echo "combined variants: $(( $(wc -l < "${GT}") - 1 )); samples: $(( $(head -1 "${GT}" | awk '{print NF}') - 2 )); tables: ${#TABLES[@]}"
    >>>

    output {
        File sample_gt_table   = "~{prefix}.sample_gt.txt"
        File sample_freq_table = "~{prefix}.sample_freq.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
