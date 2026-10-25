version 1.0

import "../../structs/Structs.wdl"

task FilterVcfForHmmIBD {

    meta {
        description: "Pre-filter a VCF for IBD inference: mask low-depth genotypes (FMT/DP < min_depth), optionally preserve the original allele-frequency annotations, recompute AN/AC/AF (optionally per population), and drop alleles/sites that are no longer represented. With output_format='bcf' (default) emits a compressed BCF ready for hmmibd-rs --from-bcf. With output_format='hmmibd_table' it FUSES the table conversion into this task: the filtered BCF is converted in-place to the hmmIBD genotype table with bcf2hmmibd (hmmibd-rs's own conversion code, byte-identical to --from-bcf; dominant-allele/first-ploidy) and the BCF is deleted, so the (often enormous) BCF is never delocalized to GCS or shipped to a separate conversion task. hmmibd-rs computes allele frequencies from the table itself, so no frequency file is produced."

        tool:          "bcftools"
        tool_version:  "1.22"
        tool_url:      "https://www.htslib.org/"
        tool_citation: "Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008."

        author: "Jonn Smith"

        outputs: {
            filtered_bcf:       "Filtered, AN/AC/AF-recomputed callset as a compressed BCF (glob: one file when output_format='bcf', empty when 'hmmibd_table')",
            filtered_bcf_index: "CSI index for filtered_bcf (same emptiness as filtered_bcf)",
            sample_gt_table:    "hmmIBD genotype table (glob: one file when output_format='hmmibd_table', else empty)"
        }
    }

    parameter_meta {
        input_vcf:             "VCF/BCF to filter prior to IBD inference. (required)"
        prefix:                "Basename for the output BCF. (required)"
        regions_bed:           "Optional BED or positions file of sites/regions to KEEP, applied FIRST — before all other filtering — via `bcftools view -T` (streamed, so no index is needed; contig names must match the VCF). When provided, the input is restricted to it; when omitted, NO region restriction is done. Use it to confine IBD to the core genome (e.g. the Pf3D7 core BED regions-20130225.Core.bed) or to a curated site list. A '.bed' filename is read as BED (0-based, half-open); any other name as 1-based CHROM<TAB>POS positions. (default: none = no restriction)"
        min_depth:             "Genotypes with FORMAT/DP below this value are set to missing (bcftools filter -S .). NOTE: this masks GT, so in dominant-allele mode (which ignores GT) it does not gate the calls — dom_min_depth is the read-depth floor there. (default: 8)"
        keep_original_af:      "If true, rename the existing INFO annotations via rename_annots_tsv before recomputing AN/AC/AF, so the original frequencies are preserved under new tags. (default: false)"
        variant_types:         "Which variant types to keep (bcftools view -v): 'snps', 'indels', or 'both'. Default 'both' keeps SNP+indel alleles (hmmibd-rs treats them on equal footing). NOTE: with split_multiallelics=false a mixed site keeps its indel alleles; strict SNP-only-ALLELES needs split_multiallelics=true. (default: both)"
        biallelic_only:        "Restrict to biallelic sites (bcftools view -m2 -M2). Default false so MULTIALLELIC sites are kept (hmmibd-rs handles up to --max-all alleles); set true (with split_multiallelics=true) for the strict biallelic path. (default: false)"
        split_multiallelics:   "Split multiallelic records into biallelic ones (bcftools norm -m-any). Default FALSE: skip the (single-threaded, slow) split and feed multiallelic sites directly to hmmibd-rs. Set true only for the strict biallelic path (needed to cleanly separate SNP alleles from indel alleles at mixed sites). (default: false)"
        max_variants:          "Optional cap on the number of variants kept. When >0 and the filtered callset exceeds it, sites are thinned evenly across the genome down to at most this many (not truncated to the first N). (default: 0 = no limit)"
        max_alt:               "Optional pre-norm site width cap. When >0, sites whose ALT count (after trimming unobserved alleles) still exceeds this are dropped BEFORE `bcftools norm` splits them — this bounds norm's per-record memory on Pf hyper-multiallelic (var-gene / indel) sites that otherwise OOM-kill it on whole-cohort joint calls. Lossy at the site level (drops those sites' SNPs); genuine multiallelic SNP sites have <=3 ALTs, so ~6-8 is a safe cap for Pf. (default: 0 = no cap)"
        mask_gq0_genotypes:    "Also set GQ0 genotypes to missing (in addition to the FORMAT/DP < min_depth mask). Older GATK (pre-4.6.0.0 GenotypeGVCFs; also GnarlyGenotyper) emitted no-/low-confidence hom-refs as 0/0 with GQ=0 instead of ./. (GATK issue #7792, fixed in PR #8741). Turning this on reproduces that fix downstream — the DP mask alone misses GQ0 calls whose DP >= min_depth. NOTE: this is the conservative choice — it drops ALL GQ0 hom-refs, including any that are genuinely well-covered ref; use a GVCF cross-reference if you need to keep those. This masks GT only, so it feeds the GT-based path (FilterVcfForHmmIBD's first-ploidy table, plus AC/AF and site/allele selection); the dominant-allele table reads FORMAT/AD and ignores GT, so it is a no-op for those call values there — the dom_min_depth/dom_min_ratio/dom_min_r1_r2 AD gates handle GQ0 hom-refs evidence-based (matching hmmibd-rs read_dom). (default: false)"
        populations_file:      "Optional sample-to-population file passed to `bcftools +fill-tags -S`; when given, AN/AC/AF are computed per population. When omitted, tags are computed across all samples. (default: none)"
        rename_annots_tsv:     "Required when keep_original_af=true: two-column TSV of old-name<TAB>new-name passed to `bcftools annotate --rename-annots`. (default: none)"
        extra_args:            "Additional `bcftools view` filters appended to the FINAL view (after AN/AC/AF are recomputed, so INFO-based expressions like MAF work). Whitespace-separated tokens; write expressions without internal spaces, e.g. \"-e MAF<0.01\" or \"-i F_MISSING<0.1\". Note the final view already applies -i 'MAX(INFO/AC) > 0'; supply extra exclusions with -e (a second -i would override the built-in one). (default: empty)"
        output_format:         "'bcf' (default) delocalizes the filtered BCF; 'hmmibd_table' converts it in-task to the hmmIBD genotype table + freq file and deletes the BCF (never delocalized). (default: bcf)"
        gt_mode:               "hmmibd_table only: 'dominant-allele' (default, from FORMAT/AD, reproduces hmmibd-rs read_dom) or 'first-ploidy' (GT[0]). (default: dominant-allele)"
        dom_min_depth:         "hmmibd_table dominant-allele: accept a call only if total AD depth > this (7 => >=8 reads). (default: 7)"
        dom_min_ratio:         "hmmibd_table dominant-allele: accept only if major/total >= this. (default: 0.7)"
        dom_min_r1_r2:         "hmmibd_table dominant-allele: accept only if minor/major < 1/this. (default: 3.0)"
        min_maf:               "hmmibd_table site prune: drop sites with minor-allele freq (of final calls) below this; >0 also drops monomorphic sites. (default: 0.01)"
        min_site_nonmissing:   "hmmibd_table site prune: drop sites with non-missing-call fraction below this. (default: 0.3)"
        max_all:               "hmmibd_table only: sites with more than this many alleles are DROPPED during conversion (whole site, never truncated) — must match hmmibd-rs --max-all so the table never exceeds what it can read (i8 call storage, hard ceiling 127). (default: 64)"
        keep_filtered_bcf:     "hmmibd_table only: also keep and delocalize the filtered BCF (+ index) instead of deleting it after the table is built. (default: false)"
        emit_progress:         "Write progress bars (percent, rate, elapsed, ETA) to stderr via pv: a byte-based bar over the input during filtration, and a windows-done bar over the parallel table conversion (output_format=hmmibd_table). pv ships in lr-basic:0.1.4 (installed at runtime if missing). Off by default. (default: false)"
        runtime_attr_override: "Override the default runtime attributes. (default: none)"
    }

    input {
        File input_vcf
        String prefix

        # Recommended default (Pf3D7 core genome):
        #   gs://broad-malaria-public/short_read_workspace_data/regions/regions-20130225.Core.bed
        File? regions_bed

        Int min_depth = 8
        Boolean keep_original_af = false
        String variant_types = "both"
        Boolean biallelic_only = false
        Boolean split_multiallelics = false
        Int max_variants = 0
        Int max_alt = 0
        Boolean mask_gq0_genotypes = false

        File? populations_file
        File? rename_annots_tsv

        String extra_args = ""
        Boolean emit_progress = false

        String output_format = "bcf"
        String gt_mode = "dominant-allele"
        Int dom_min_depth = 7
        Float dom_min_ratio = 0.7
        Float dom_min_r1_r2 = 3.0
        Float min_maf = 0.01
        Float min_site_nonmissing = 0.3
        Int max_all = 64
        Boolean keep_filtered_bcf = false

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

        # Pre-norm shrink to keep `norm` memory bounded on wide multiallelic sites (the Pf var-gene /
        # indel blowups that OOM-kill norm on whole-cohort joint calls). `norm -m-any` memory scales
        # with (samples x alleles) because it splits a site into one biallelic record per ALT while
        # carrying FORMAT/AD (Number=R). Two cheap streaming stages cut the allele width BEFORE norm:
        #   1. --trim-alt-alleles: drop ALT alleles not seen in any (DP-masked) genotype. Lossless --
        #      those alleles have AC=0 and are trimmed at the end anyway, so the output is identical;
        #      only the intermediate (and thus norm's per-record work) shrinks.
        #   2. optional N_ALT cap: drop sites whose remaining ALT count still exceeds max_alt, so a
        #      pathological hyper-multiallelic site never reaches norm's per-allele split at all.
        # Note: memory here is per-RECORD (widest site), so this width cut -- not genomic sharding --
        # is what bounds norm's peak memory.
        TRIM_STEP=(bcftools view --trim-alt-alleles -Ou)
        MAXALT_STEP=(cat)
        if [[ ~{max_alt} -gt 0 ]] ; then MAXALT_STEP=(bcftools view -e "N_ALT > ~{max_alt}" -Ou) ; fi

        # Optional region/site restriction, applied FIRST (before every other stage), so unwanted
        # regions -- subtelomeric / hypervariable / VAR-gene, or anything outside a curated site
        # list -- never enter IBD. Applied iff regions_bed is provided. Streamed with
        # `bcftools view -T` (targets file), which filters the piped input by position without an
        # index (unlike -R, which needs random access). `cat` is a no-op passthrough when no
        # regions_bed is given.
        REGIONS_STEP=(cat)
        if [[ -n "~{regions_bed}" ]] ; then
            REGIONS_STEP=(bcftools view -T "~{regions_bed}" -Ou)
        fi

        # Optional progress bar: stream the input through pv so its byte-based percent, rate, elapsed,
        # and ETA go to stderr. The pipe applies backpressure, so pv's read rate tracks the whole
        # filter's end-to-end throughput (not just decompression). pv is installed if missing (best
        # effort; the run continues without a bar on failure). `cat` is a no-op passthrough when off.
        PV=(cat)
        if [[ "~{emit_progress}" == "true" ]]; then
            command -v pv >/dev/null 2>&1 || { { apt-get update -qq && apt-get install -y -qq pv ; } >/dev/null 2>&1 || echo "WARN: pv unavailable and install failed; continuing without a progress bar" >&2 ; }
            if command -v pv >/dev/null 2>&1 ; then PV=(pv -f -p -t -e -r -b -i 15 -N filter) ; fi
        fi

        # keeping the original allele frequencies requires a rename map
        if [[ "~{keep_original_af}" == "true" && -z "~{rename_annots_tsv}" ]] ; then
            echo "ERROR: keep_original_af=true requires rename_annots_tsv to be provided." >&2
            exit 1
        fi

        # Mask low-depth genotypes, (optionally) preserve original AF via rename, recompute
        # AN/AC/AF (optionally per population), then trim now-unrepresented alt alleles/sites.
        if [[ "~{keep_original_af}" == "true" ]] ; then
            "${PV[@]}" ~{input_vcf} \
              | "${REGIONS_STEP[@]}" \
              | bcftools annotate -x '~{fmt_keep}' -Ou - \
              | bcftools filter -S . -e "~{gt_mask_expr}" -Ou \
              | bcftools view "${TYPE_FLAGS[@]}" -Ou \
              | "${TRIM_STEP[@]}" \
              | "${MAXALT_STEP[@]}" \
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
            "${PV[@]}" ~{input_vcf} \
              | "${REGIONS_STEP[@]}" \
              | bcftools annotate -x 'INFO,~{fmt_keep}' -Ou - \
              | bcftools filter -S . -e "~{gt_mask_expr}" -Ou \
              | bcftools view "${TYPE_FLAGS[@]}" -Ou \
              | "${TRIM_STEP[@]}" \
              | "${MAXALT_STEP[@]}" \
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
        # Table generation (output_format=hmmibd_table): convert the LOCAL filtered BCF to the hmmIBD
        # genotype table with bcf2hmmibd and DELETE the BCF, so the enormous file is never delocalized
        # or shipped to a second task. bcf2hmmibd is hmmibd-rs's OWN conversion code
        # (BcfGenotype::new_from_processing_bcf -> into_genotype_siteinfo -- the exact path
        # `hmmibd-rs --from-bcf` runs), so the per-shard table is byte-identical to what --from-bcf
        # would build from this BCF under the same filter config. See
        # docker/lr-hmmibd-rs/bcf2hmmibd/README.md for the equivalence proof.
        if [[ "~{output_format}" == "hmmibd_table" ]]; then
            GT="~{prefix}.sample_gt.txt"
            BCF="~{prefix}.filtered.bcf"

            # Allele-width cap (matches the old conversion's nA>max_all drop): hmmibd-rs stores each
            # call as i8 and sizes its per-site frequency matrix by --max-all, so a site with more than
            # max_all alleles must not reach it. Drop those sites here (N_ALT = nA-1, so nA>max_all means
            # N_ALT>=max_all). This is a filtration choice on the BCF, so bcf2hmmibd stays byte-identical
            # to `hmmibd-rs --from-bcf` run on the SAME (capped) BCF.
            if [[ ~{max_all} -gt 0 ]]; then
                bcftools view -e "N_ALT >= ~{max_all}" "${BCF}" -Ob --threads "${NUM_CPUS}" -o "${BCF}.capped"
                mv "${BCF}.capped" "${BCF}"
            fi

            # Map the WDL filter knobs to a BcfFilterArgs TOML (read_dom AD gates + per-site filters).
            # nonmissing_rates is deliberately EMPTY: that iterative site+SAMPLE missingness refinement
            # needs whole-callset missingness, which is unavailable when the callset is sharded across
            # input VCFs (a sample dropped in one shard but not another would desync the per-shard
            # sample columns and break the stitch). The per-SITE filters (PASS, nallele>1,
            # min_site_nonmissing, min_maf) are computed over ALL samples present in each shard, so they
            # ARE shard-consistent and stay on. filter_column_pass_only=true applies the PASS-only site
            # filter that --from-bcf applies (the old awk conversion omitted it).
            cat > bcf_filter.toml <<TOML
min_depth = ~{dom_min_depth}
min_ratio = ~{dom_min_ratio}
min_r1_r2 = ~{dom_min_r1_r2}
min_maf = ~{min_maf}
min_site_nonmissing = ~{min_site_nonmissing}
filter_column_pass_only = true
major_minor_alleles_must_be_snps = true
nonmissing_rates = []
TOML

            # --from-bcf is a required no-op flag here (hmmibd-rs's shared CLI demands it whenever
            # --bcf-filter-config is given; bcf2hmmibd always converts regardless of it). --min-snp-sep 1
            # keeps the per-shard table maximally dense (drops only exact-duplicate positions); the
            # MAF-prioritized spacing/density thinning runs once on the COMBINED callset in
            # StitchHmmIBDTables, and hmmibd-rs re-applies --min-snp-sep at IBD time. chrom is emitted as
            # the BCF contig NAME (hmmibd-rs reads string chrom; the stitch version-sorts it).
            bcf2hmmibd \
                --from-bcf \
                -i "${BCF}" \
                -o "${GT}" \
                --bcf-read-mode "~{gt_mode}" \
                --bcf-filter-config bcf_filter.toml \
                --min-snp-sep 1
            rm -f bcf_filter.toml

            # Delete the BCF (never delocalized) unless the caller asked to keep it as an output.
            if [[ "~{keep_filtered_bcf}" != "true" ]]; then rm -f ~{prefix}.filtered.bcf ~{prefix}.filtered.bcf.csi ; fi
            echo "hmmIBD table (bcf2hmmibd): variants=$(( $(wc -l < "${GT}") - 1 )); samples=$(( $(head -1 "${GT}" | awk '{print NF}') - 2 )); mode=~{gt_mode}; keep_bcf=~{keep_filtered_bcf}"
        fi
    >>>

    output {
        # glob-backed so each is present only for the chosen output_format: 'bcf' delocalizes the BCF
        # (+ index); 'hmmibd_table' deletes the BCF and delocalizes only the small genotype table.
        Array[File] filtered_bcf       = glob("~{prefix}.filtered.bcf")
        Array[File] filtered_bcf_index = glob("~{prefix}.filtered.bcf.csi")
        Array[File] sample_gt_table    = glob("~{prefix}.sample_gt.txt")
    }

    #########################
    # mem_gb default 16: `bcftools norm -m-any` (split multiallelics) is the memory hog on
    # whole-cohort joint calls — Pf sites can carry hundreds of ALT alleles and FORMAT/AD is
    # Number=R, so per-record memory scales with (samples x alleles). A 4 GB default OOM-kills
    # norm on realistic malaria joint calls (e.g. a 28 GB single-contig cohort VCF). Override
    # runtime_attr_override.mem_gb higher (32-64) for the widest cohorts; NOTE Cromwell's
    # memory-retry only bumps memory across attempts if the workspace sets memory_retry_multiplier.
    # preemptible_tries 0: this task is long (hours) on whole-cohort joint calls, and a Spot/preemptible
    # VM eviction restarts it from scratch — observed ~3 h shards preempted repeatedly and never
    # converging. On-demand costs more per hour but actually finishes. Override to >0 only for inputs
    # small enough to complete inside a preemption window.
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-hmmibd-rs:0.1.8"
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

        max_all:               "Maximum number of unique alleles per site (--max-all). Allele calls are stored as i8, so the hard ceiling is 127; memory scales linearly (per-site freq matrix). (default: 64)"
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
        ibd_output_disk_factor: "Disk multiplier applied to the UNCOMPRESSED input genotype volume to provision for the O(N^2) hmm.txt output (disk_size = 50 + ibd_output_disk_factor*genotype_vol_gb + 3*input). Large by design (hmm.txt was ~175x the table on a 32k-sample cohort); raise for very large N, or set runtime_attr_override.disk_gb directly. (default: 200)"
        bcf_decompression_factor: "When from_bcf=true the input is a compressed BCF; multiply its file size by this to estimate the uncompressed genotype volume the output scales with. Ignored when from_bcf=false (text table already uncompressed). (default: 8.0)"
        ibd_core_seconds_per_gb2: "CPU auto-sizing calibration: hmmibd-rs core-seconds per (uncompressed genotype GB)^2 for the O(N^2) pair sweep. From a Gambia run (8 cores, 2.68 GB table, ~18.9 days). cpu_cores is derived so the run finishes within ibd_target_runtime_hours. (default: 1.82e6)"
        ibd_target_runtime_hours: "CPU auto-sizing target wall-clock (hours). cpu_cores = smallest power of two (<= ibd_max_cores) that finishes the estimated work within this time. (default: 48)"
        ibd_max_cores: "Cap on auto-sized cpu_cores (power of two). Above ~3.5 GB genotype volume even this cannot hit the target (O(N^2) output) -> pair-shard. (default: 128)"

        extra_args:            "Additional command-line args appended verbatim to the hmmibd-rs invocation. (default: empty)"
        emit_progress:         "Pass hmmibd-rs --print-progress so it periodically reports progress during the O(N^2) pair sweep. Each PROGRESS line is reformatted on stderr to prepend a wall-clock UTC timestamp and append an estimated-time-remaining column, e.g. `[14:03:27 UTC] PROGRESS 43.7%  pairs=1830412  elapsed=00:10:12  ETA=00:13:08` (ETA = elapsed*(1-frac)/frac, linear extrapolation). (default: false)"
        runtime_attr_override: "Override the default runtime attributes. (default: none)"
    }

    input {
        File input_bcf
        String prefix

        File? data_file2
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

        Int max_all = 64
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

        # Disk sizing for the O(N^2) hmm.txt output (see disk_size below).
        Int ibd_output_disk_factor = 200
        Float bcf_decompression_factor = 8.0

        # CPU auto-sizing for the O(N^2) IBD sweep (see cpu_cores below): enough cores to finish
        # within ibd_target_runtime_hours, rounded UP to a power of two and capped at ibd_max_cores.
        Float ibd_core_seconds_per_gb2 = 1820000.0
        Int ibd_target_runtime_hours = 48
        Int ibd_max_cores = 128

        String extra_args = ""
        Boolean emit_progress = false

        RuntimeAttr? runtime_attr_override
    }

    # hmmibd-rs OUTPUT scales with sample-PAIRS (~N^2/2), NOT with the genotype-table input:
    # .hmm.txt is >=1 UNCOMPRESSED segment line per pair and is the disk hog, and .hmm_fract.txt is
    # ~1 line/pair. A StorageFull panic here -- hours/days into the O(N^2) sweep -- is the usual
    # failure (a ~32k-sample Gambia run filled a 104 GB disk at ~21% of pairs; it needed ~470 GB).
    # We can't read N from the file at workflow-graph time, so size disk off the input GENOTYPE
    # VOLUME, compression-aware, times a generous factor for the quadratic output:
    #   from_bcf=false -> input is the UNCOMPRESSED hmmIBD text table: use its size as-is
    #   from_bcf=true  -> input is a COMPRESSED BCF: scale up by bcf_decompression_factor first
    # ibd_output_disk_factor is large on purpose: hmm.txt was ~175x the uncompressed table on that
    # 32k-sample cohort, and the ratio GROWS with N (it is O(N^2) output over O(N) input). For very
    # large N raise ibd_output_disk_factor (or set runtime_attr_override.disk_gb), or cut the output
    # with filt_ibd_only / filt_min_seg_cm / suppress_frac.
    Float input_data_gb   = size(input_bcf, "GB") + size(data_file2, "GB")
    Float genotype_vol_gb = if from_bcf then input_data_gb * bcf_decompression_factor else input_data_gb
    Int disk_size = 50 + ceil(ibd_output_disk_factor * genotype_vol_gb) + ceil(3.0 * input_data_gb)

    # CPU auto-sizing. hmmibd-rs runtime is dominated by the O(N^2) pair sweep, cost ~ (pairs x sites)
    # = O(N^2 x S). The genotype table is ~ N x S bytes and the pipeline thins S to a bounded
    # genome-wide density, so total core-seconds scale ~ with the (uncompressed) table size SQUARED.
    # Calibrated on a Gambia run: 8 cores, 2.68 GB table, ~526M pairs, ~18.9 days => ~13.1M
    # core-seconds => ~1.82e6 core-sec / GB^2. cores_needed = core_sec / target_seconds; round UP to a
    # power of two, clamp to [8, ibd_max_cores]. NOTE: output is O(N^2), so above ~3.5 GB even 128
    # cores cannot hit 2 days -> pair-shard the run, or raise ibd_target_runtime_hours / ibd_max_cores.
    Float ibd_core_seconds = ibd_core_seconds_per_gb2 * genotype_vol_gb * genotype_vol_gb
    Float ibd_raw_cores    = ibd_core_seconds / (ibd_target_runtime_hours * 3600)
    Int ibd_pow2_cores = if      ibd_raw_cores <= 8.0  then 8
                         else if ibd_raw_cores <= 16.0 then 16
                         else if ibd_raw_cores <= 32.0 then 32
                         else if ibd_raw_cores <= 64.0 then 64
                         else 128
    Int ibd_cpu_cores = if ibd_pow2_cores < ibd_max_cores then ibd_pow2_cores else ibd_max_cores
    # keep the GCP custom-machine mem/CPU ratio valid (>= 0.9 GB/vCPU): ~1 GB/core, floor 16.
    Int ibd_mem_gb = if ibd_cpu_cores < 16 then 16 else ibd_cpu_cores

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

        # Reformat hmmibd-rs --print-progress lines on stderr: prepend a wall-clock (UTC) timestamp
        # and append an estimated-time-remaining column, derived from the fraction + elapsed seconds
        # the binary already prints (ETA = elapsed * (1 - frac) / frac). Pure bash + coreutils `date`
        # -- the hmmibd-rs image is debian-slim with no awk. Non-PROGRESS stderr (the Arguments dump,
        # any error trace) passes through verbatim. `set +x` keeps the per-line loop out of xtrace.
        _fmt_hms() { local s=$1 h m; if [ "${s}" -lt 0 ] 2>/dev/null; then printf '?'; return; fi; h=$((s/3600)); m=$(((s%3600)/60)); s=$((s%60)); printf '%02d:%02d:%02d' "${h}" "${m}" "${s}"; }
        _hmm_progress() {
            set +x
            local line frac pairs el ip dp pm eta
            while IFS= read -r line; do
                case "${line}" in
                    PROGRESS*)
                        read -r frac pairs el _ <<< "${line#PROGRESS}"   # frac pairs elapsed_s
                        ip=${frac%%.*}; dp=${frac#*.}
                        if [ "${dp}" = "${frac}" ]; then dp=0; fi
                        dp=${dp}000; dp=${dp:0:3}
                        pm=$(( 10#${ip:-0} * 1000 + 10#${dp:-0} ))        # progress in permille (0..1000)
                        if [ "${pm}" -gt 0 ]; then eta=$(( el * (1000 - pm) / pm )); else eta=-1; fi
                        printf '[%s UTC] PROGRESS %s.%s%%  pairs=%s  elapsed=%s  ETA=%s\n' \
                            "$(date -u +%H:%M:%S)" "$(( pm / 10 ))" "$(( pm % 10 ))" \
                            "${pairs}" "$(_fmt_hms "${el}")" "$(_fmt_hms "${eta}")"
                        ;;
                    *) printf '%s\n' "${line}" ;;
                esac
            done
        }

        hmmibd-rs \
            ~{true="--from-bcf" false="" from_bcf} \
            -i ~{input_bcf} \
            -o ~{prefix} \
            ~{read_mode_arg} \
            ~{"-I " + data_file2} \
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
            ~{true="--print-progress" false="" emit_progress} \
            ~{true="--suppress-frac" false="" suppress_frac} \
            ~{true="--bcf-to-bin-file" false="" bcf_to_bin_file} \
            ~{true="--bcf-to-bin-file-by-chromosome" false="" bcf_to_bin_file_by_chromosome} \
            ~{extra_args} \
            2> >(_hmm_progress >&2)
    >>>

    output {
        # Glob so each output exists only when the chosen flags actually write it.
        Array[File] ibd_segments     = glob("~{prefix}.hmm.txt")
        Array[File] ibd_fraction     = glob("~{prefix}.hmm_fract.txt")
        Array[File] binary_genotypes = glob("~{prefix}*.bin")
    }

    #########################
    # IBD is O(N^2) over sample PAIRS (e.g. ~12k samples -> ~69M pairs), and hmmibd-rs parallelizes
    # over pairs, so cpu_cores is the main throughput lever. cpu_cores/mem_gb are AUTO-SIZED above
    # (ibd_cpu_cores/ibd_mem_gb) to finish within ibd_target_runtime_hours (default 48h); override via
    # runtime_attr_override.cpu_cores (+ mem_gb) to force a size, or consider par_mode=1 / pair-sharding
    # for very large sample sets. preemptible_tries=0 because this can run for days; a Spot eviction
    # wastes the whole attempt.
    RuntimeAttr default_attr = object {
        cpu_cores:          ibd_cpu_cores,
        mem_gb:             ibd_mem_gb,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-hmmibd-rs:0.1.8"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.4"
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
        description: "NOTE: standalone BIALLELIC-only converter, no longer used by the pipelines (FilterVcfForHmmIBD fuses a multiallelic-aware conversion into its task). Kept for standalone/biallelic use; for multiallelic use the fused FilterVcfForHmmIBD output_format='hmmibd_table' path. Convert a bi-allelic BCF/VCF to the hmmIBD text genotype table: a tab-delimited matrix with columns chrom(int) pos then one call per sample, alleles coded 0/1 and -1 for missing. The per-sample call is derived per gt_mode: 'dominant-allele' (default) reproduces hmmibd-rs --bcf-read-mode dominant-allele (src/bcf.rs read_dom) — from FORMAT/AD it takes the highest-depth allele (ties -> lower allele index) and keeps it only when total>dom_min_depth & major/total>=dom_min_ratio & minor/major<1/dom_min_r1_r2, else -1; GT is intentionally ignored, matching hmmibd-rs, so this is the right choice for polyclonal Pf (majority-clone allele, not GATK's ref-biased GT[0]). 'first-ploidy' instead takes the first allele of GT. Because it runs on the filtered bi-allelic BCF, dominant-allele here matches an hmmibd-rs --from-bcf dominant-allele run on that same BCF. Sites are then pruned on the FINAL calls to match hmmibd-rs --from-bcf site filtering (min_maf, min_site_nonmissing), so the table and --from-bcf paths use the same sites; min_maf>0 also drops monomorphic / no-alt-expressed sites, which carry no IBD information. hmmibd-rs computes allele frequencies from the table itself, so no frequency file is produced. Non-numeric contigs (MIT/API) are dropped because hmmIBD requires integer chromosomes."

        tool:          "bcftools"
        tool_version:  "1.22"
        tool_url:      "https://www.htslib.org/"

        author: "Jonn Smith"

        outputs: {
            sample_gt_table:   "hmmIBD text genotype table (<prefix>.sample_gt.txt); hmmibd-rs input in text mode (from_bcf=false)"
        }
    }

    parameter_meta {
        input_bcf:             "Bi-allelic, SNP-filtered BCF/VCF (as produced by FilterVcfForHmmIBD) to convert. Must carry FORMAT/AD when gt_mode='dominant-allele'. (required)"
        prefix:                "Basename for the output tables. (required)"
        gt_mode:               "How to derive each per-sample call: 'dominant-allele' (default) = max-depth allele from FORMAT/AD with hmmibd-rs read_dom gating (the right choice for polyclonal Pf); 'first-ploidy' = first allele of GT. (default: dominant-allele)"
        dom_min_depth:         "dominant-allele only: a call is accepted only if total AD depth > dom_min_depth, else missing. The gate is strict '>', so 7 requires >=8 reads. (default: 7; diverges from hmmibd-rs's default of 5)"
        dom_min_ratio:         "dominant-allele only: minimum major-allele fraction (major/total >= dom_min_ratio) to accept a call. Matches hmmibd-rs min_ratio. (default: 0.7)"
        dom_min_r1_r2:         "dominant-allele only: a call is accepted only if minor/major < 1/dom_min_r1_r2 (i.e. the minor allele is well below the major). Matches hmmibd-rs min_r1_r2. (default: 3.0)"
        min_maf:               "Drop sites whose minor-allele frequency among non-missing FINAL calls is below this (matches hmmibd-rs min_maf). Values >0 also drop monomorphic / no-alt-expressed sites (no IBD information); set 0 to keep rare-variant sites (which carry strong IBD signal). (default: 0.01)"
        min_site_nonmissing:   "Drop sites where the fraction of samples with a non-missing call is below this (matches hmmibd-rs min_site_nonmissing); sparse sites give unreliable allele frequencies. (default: 0.3)"
        runtime_attr_override: "Override the default runtime attributes. (default: none)"
    }

    input {
        File input_bcf
        String prefix

        String gt_mode = "dominant-allele"
        Int dom_min_depth = 7
        Float dom_min_ratio = 0.7
        Float dom_min_r1_r2 = 3.0
        Float min_maf = 0.01
        Float min_site_nonmissing = 0.3

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(input_bcf, "GB"))

    command <<<
        set -euxo pipefail

        GT="~{prefix}.sample_gt.txt"
        MODE="~{gt_mode}"

        # header: chrom  pos  <sample1> <sample2> ...
        { printf 'chrom\tpos'; bcftools query -l "~{input_bcf}" | awk '{printf "\t%s", $0}'; printf '\n'; } > "${GT}"

        # body: integer chrom, pos, one call per sample (-1 missing). Contig prefix stripped to int.
        case "${MODE}" in
            dominant-allele)
                # Reproduce hmmibd-rs --bcf-read-mode dominant-allele (src/bcf.rs read_dom): from the
                # bi-allelic FORMAT/AD (ref,alt), pick the higher-depth allele (tie -> lower index),
                # and accept it only when total>dom_min_depth & major/total>=dom_min_ratio &
                # minor/major < 1/dom_min_r1_r2; otherwise the site is missing (-1). GT is ignored
                # (matches hmmibd-rs). Same math as read_dom on the same bi-allelic BCF.
                bcftools query -f '%CHROM\t%POS[\t%AD]\n' "~{input_bcf}" \
                | awk -v md=~{dom_min_depth} -v mr=~{dom_min_ratio} -v mrr=~{dom_min_r1_r2} 'BEGIN{FS=OFS="\t"; inv=1.0/mrr}
                {
                    c=$1; sub(/^Pf3D7_/,"",c); sub(/_v[0-9]+$/,"",c)
                    if (c !~ /^[0-9]+$/) next
                    printf "%d\t%d", c+0, $2
                    for (i=3;i<=NF;i++){
                        g=-1
                        if ($i!="." && $i!=""){
                            split($i,a,","); ref=(a[1]=="."?0:a[1]+0); alt=(a[2]=="."?0:a[2]+0)
                            total=ref+alt
                            if (alt>ref){ major=alt; minor=ref; ma=1 } else { major=ref; minor=alt; ma=0 }
                            if (total>md && major/total>=mr && minor/major<inv) g=ma
                        }
                        printf "\t%s", g
                    }
                    printf "\n"
                }' > body_raw.txt
                ;;
            first-ploidy)
                # First allele of GT per sample (-1 missing); matches hmmibd-rs --bcf-read-mode first-ploidy.
                bcftools query -f '%CHROM\t%POS[\t%GT]\n' "~{input_bcf}" | awk 'BEGIN{FS=OFS="\t"}
                {
                    c=$1; sub(/^Pf3D7_/,"",c); sub(/_v[0-9]+$/,"",c)
                    if (c !~ /^[0-9]+$/) next
                    printf "%d\t%d", c+0, $2
                    for (i=3;i<=NF;i++){ n=split($i,a,/[\/|]/); g=a[1]; if (g=="."||g=="") g=-1; printf "\t%s", g }
                    printf "\n"
                }' > body_raw.txt
                ;;
            *)
                echo "ERROR: gt_mode must be 'dominant-allele' or 'first-ploidy' (each-ploidy is only available via hmmibd-rs --from-bcf)." >&2
                exit 1
                ;;
        esac

        # Site-level prune (match hmmibd-rs --from-bcf so the table and --from-bcf paths agree on
        # which sites are used): drop sites below min_site_nonmissing (fraction of samples with a
        # non-missing call -- sparse sites give unreliable allele frequencies) and below min_maf
        # (minor-allele frequency among non-missing calls; min_maf>0 also drops monomorphic /
        # no-alt-expressed sites, which carry no IBD information). Computed on the FINAL calls, so it
        # honors dominant-allele, over the full sample set (each site is present here with all samples).
        awk -v mm=~{min_maf} -v mn=~{min_site_nonmissing} 'BEGIN{FS=OFS="\t"}
        {
            ns=NF-2; c0=0; c1=0; nm=0
            for (i=3;i<=NF;i++){ if($i=="0"){c0++;nm++} else if($i=="1"){c1++;nm++} }
            if (ns<=0 || nm==0) next                       # no samples / all missing -> drop
            if ((nm/ns) < mn) next                         # too many missing calls
            maf=(c0<c1?c0:c1)/nm
            if (maf < mm) next                             # below MAF (min_maf>0 also drops monomorphic)
            print
        }' body_raw.txt >> "${GT}"
        rm -f body_raw.txt

        echo "variants: $(( $(wc -l < "${GT}") - 1 )); samples: $(( $(head -1 "${GT}" | awk '{print NF}') - 2 )); mode: ${MODE}"
    >>>

    output {
        File sample_gt_table   = "~{prefix}.sample_gt.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.4"
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
        description: "Combine many per-region hmmIBD genotype tables (one per input VCF, as produced by BcfToSampleTable) into a SINGLE hmmIBD genotype table for hmmibd-rs. This is the scalable combine step: each huge input VCF is first reduced to its compact integer genotype table, and only these small tables are stitched here — a merged multi-hundred-GB BCF is never materialized. Assumes every input table carries the SAME samples in the SAME column order (i.e. the input VCFs are disjoint regions of ONE joint call, split across files); the sample headers are verified identical and the task aborts if they differ. Rows (sites) from all tables are concatenated and coordinate-sorted by chrom (version-sorted, so integer or contig-name chrom both order correctly) then pos (hmmibd-rs assumes sorted input); closely spaced / linked variants are then thinned (keep at most thin_max_per_window sites per thin_window_bp window and enforce thin_min_snp_sep bp minimum spacing, prioritizing higher minor-allele frequency) because hmmIBD assumes markers are ~independent given IBD state and over-dense linked SNPs inflate IBD. hmmibd-rs computes allele frequencies from the combined table itself, so no frequency file is produced. An optional max_variants cap is applied after thinning to the combined, genome-wide callset (even stride, not truncated). Input regions are expected to be non-overlapping."

        tool:          "bcftools"
        tool_version:  "1.22"
        tool_url:      "https://www.htslib.org/"

        author: "Jonn Smith"

        outputs: {
            sample_gt_table:   "Combined hmmIBD text genotype table (<prefix>.sample_gt.txt): all input sites, coordinate-sorted, one first-allele call per sample (-1 missing)"
        }
    }

    parameter_meta {
        gt_tables:             "hmmIBD genotype tables to combine, one per input VCF (from BcfToSampleTable / FilterVcfForHmmIBD output_format='hmmibd_table'). All must share the same samples in the same column order. (required)"
        prefix:                "Basename for the combined output tables. (required)"
        thin_window_bp:        "LD/density thinning window size in bp (SLIDING window, applied by bcftools +prune). Keep at most thin_max_per_window sites per window. (default: 2000)"
        thin_max_per_window:   "LD/density thinning: max sites kept per thin_window_bp window, chosen by highest MAF. hmmIBD assumes markers ~independent given IBD state; over-dense linked SNPs inflate IBD. 0 disables the density cap. (default: 12)"
        thin_min_snp_sep:      "LD/density thinning: minimum bp spacing between kept sites (MAF-prioritized). 0 disables the spacing rule. Set both this and thin_max_per_window to 0 to skip thinning entirely. (default: 50)"
        max_variants:          "Optional cap on the number of variants in the COMBINED callset, applied AFTER thinning. When >0 and the count still exceeds it, sites are thinned evenly across the genome down to at most this many (not truncated to the first N). Applied genome-wide here (not per input file). (default: 0 = no limit)"
        runtime_attr_override: "Override the default runtime attributes. (default: none)"
    }

    input {
        Array[File] gt_tables
        String prefix
        Int thin_window_bp = 2000
        Int thin_max_per_window = 12
        Int thin_min_snp_sep = 50
        Int max_variants = 0
        RuntimeAttr? runtime_attr_override
    }

    # Sort is external (disk-backed); size for input + sort temp + output.
    Int disk_size = 20 + ceil(10.0 * size(gt_tables, "GB"))

    command <<<
        set -euxo pipefail

        GT="~{prefix}.sample_gt.txt"

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
            | sort -T . -k1,1V -k2,2n > body_sorted.txt

        # Thin closely spaced / linked variants (default on): enforce a minimum thin_min_snp_sep bp
        # spacing and keep at most thin_max_per_window sites per SLIDING thin_window_bp window,
        # PRIORITIZING higher minor-allele frequency (the more informative markers). hmmIBD assumes
        # markers are ~independent given IBD state; over-dense linked SNPs violate that and inflate IBD
        # sharing. Runs on the combined, coordinate-sorted callset (genome-wide; a window never crosses
        # a shard boundary here) using the dominant-allele-call MAF. The sliding-window density cap is
        # applied by bcftools +prune on a sites-only VCF carrying INFO/MAF, i.e. it matches a
        # `+prune -w <bp> -n <N> --AF-tag MAF -N maxAF` run; the minimum spacing is a preceding
        # MAF-greedy pass (+prune has no min-distance option). Set thin_max_per_window=0 and
        # thin_min_snp_sep=0 to disable.
        if [[ ~{thin_max_per_window} -gt 0 || ~{thin_min_snp_sep} -gt 0 ]]; then
            BEFORE=$(wc -l < body_sorted.txt)

            # (chrom, pos, MAF) for every site, from the combined calls. Multiallelic-aware:
            # MAF = 1 - major-allele frequency over all called alleles (== min(f0,f1) for biallelic).
            awk 'BEGIN{FS=OFS="\t"} { delete cnt; nm=0; mx=0
                 for(i=3;i<=NF;i++){ v=$i; if(v!="-1" && v!="" && v!="."){ cnt[v]++; nm++; if(cnt[v]>mx)mx=cnt[v] } }
                 maf=(nm==0?0:1-mx/nm); printf "%s\t%s\t%.6f\n",$1,$2,maf }' body_sorted.txt > site_maf.txt

            # Minimum spacing: greedy by MAF descending; drop a site if a higher-MAF kept site is
            # within thin_min_snp_sep bp (checked against its own +/- S-block neighbors).
            if [[ ~{thin_min_snp_sep} -gt 0 ]]; then
                sort -T . -k1,1V -k3,3gr -k2,2n site_maf.txt \
                | awk -v S=~{thin_min_snp_sep} 'BEGIN{FS=OFS="\t"}
                  { chrom=$1; pos=$2+0; reject=0; sb=int(pos/S)
                    for(d=-1;d<=1;d++){ akey=chrom SUBSEP (sb+d)
                      if(akey in acc){ m=split(acc[akey],pp," "); for(i=1;i<=m;i++){ dd=pos-pp[i]; if(dd<0)dd=-dd; if(dd<S){reject=1;break} } }
                      if(reject)break }
                    if(reject)next
                    akey=chrom SUBSEP int(pos/S); acc[akey]=acc[akey] " " pos
                    print $1,$2,$3 }' > site_maf.spaced.txt
                mv site_maf.spaced.txt site_maf.txt
            fi

            # Sliding-window density cap via bcftools +prune (keep highest-MAF per window). Build a
            # sites-only VCF (contigs from the data, INFO/MAF per site, coordinate-sorted) and prune.
            if [[ ~{thin_max_per_window} -gt 0 ]]; then
                { echo '##fileformat=VCFv4.2'
                  echo '##INFO=<ID=MAF,Number=A,Type=Float,Description="minor allele frequency">'
                  cut -f1 site_maf.txt | sort -uV | awk '{print "##contig=<ID=" $0 ">"}'
                  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                  sort -k1,1V -k2,2n site_maf.txt | awk 'BEGIN{FS=OFS="\t"}{ printf "%s\t%s\t.\tN\tA\t.\t.\tMAF=%s\n",$1,$2,$3 }'
                } > sites.vcf
                bcftools +prune -w ~{thin_window_bp}bp -n ~{thin_max_per_window} --AF-tag MAF -N maxAF sites.vcf -Ov \
                  | awk 'BEGIN{FS=OFS="\t"} !/^#/ { print $1,$2 }' > kept.txt
            else
                cut -f1,2 site_maf.txt > kept.txt
            fi

            # keep only surviving sites, preserving coordinate order
            awk 'BEGIN{FS=OFS="\t"} NR==FNR{keep[$1 SUBSEP $2]=1; next} (($1 SUBSEP $2) in keep)' kept.txt body_sorted.txt > body_thinned.txt
            mv body_thinned.txt body_sorted.txt
            rm -f site_maf.txt sites.vcf kept.txt
            echo "thinning: ${BEFORE} -> $(wc -l < body_sorted.txt) sites (<= ~{thin_max_per_window}/~{thin_window_bp}bp sliding, min ~{thin_min_snp_sep}bp, MAF-prioritized)"
        fi

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

        echo "combined variants: $(( $(wc -l < "${GT}") - 1 )); samples: $(( $(head -1 "${GT}" | awk '{print NF}') - 2 )); tables: ${#TABLES[@]}"
    >>>

    output {
        File sample_gt_table   = "~{prefix}.sample_gt.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.4"
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
