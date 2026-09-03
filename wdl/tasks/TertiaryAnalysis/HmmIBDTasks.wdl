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
        populations_file:      "Optional sample-to-population file passed to `bcftools +fill-tags -S`; when given, AN/AC/AF are computed per population. When omitted, tags are computed across all samples. (default: none)"
        rename_annots_tsv:     "Required when keep_original_af=true: two-column TSV of old-name<TAB>new-name passed to `bcftools annotate --rename-annots`. (default: none)"
        extra_args:            "Additional command-line args appended verbatim to the final bcftools view invocation. (default: empty)"
        runtime_attr_override: "Override the default runtime attributes. (default: none)"
    }

    input {
        File input_vcf
        String prefix

        Int min_depth = 5
        Boolean keep_original_af = false

        File? populations_file
        File? rename_annots_tsv

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(input_vcf, "GB"))

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

        # keeping the original allele frequencies requires a rename map
        if [[ "~{keep_original_af}" == "true" && -z "~{rename_annots_tsv}" ]] ; then
            echo "ERROR: keep_original_af=true requires rename_annots_tsv to be provided." >&2
            exit 1
        fi

        # Mask low-depth genotypes, (optionally) preserve original AF via rename, recompute
        # AN/AC/AF (optionally per population), then trim now-unrepresented alt alleles/sites.
        if [[ "~{keep_original_af}" == "true" ]] ; then
            bcftools filter -S . -e "FMT/DP < ~{min_depth}" -Ou ~{input_vcf} \
              | bcftools annotate --rename-annots ~{rename_annots_tsv} -Ou \
              | bcftools +fill-tags -Ou -- ~{"-S " + populations_file} -t AN,AC,AF \
              | bcftools view --trim-alt-alleles -i 'MAX(INFO/AC) > 0' -Ob --threads ${NUM_CPUS} ~{extra_args} -o ~{prefix}.filtered.bcf
        else
            bcftools filter -S . -e "FMT/DP < ~{min_depth}" -Ou ~{input_vcf} \
              | bcftools +fill-tags -Ou -- ~{"-S " + populations_file} -t AN,AC,AF \
              | bcftools view --trim-alt-alleles -i 'MAX(INFO/AC) > 0' -Ob --threads ${NUM_CPUS} ~{extra_args} -o ~{prefix}.filtered.bcf
        fi

        bcftools index --threads ${NUM_CPUS} ~{prefix}.filtered.bcf
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
        input_bcf:             "Filtered BCF to run IBD on (read via hmmibd-rs --from-bcf; requires FORMAT/AD under the default read mode). (required)"
        prefix:                "Output prefix; produces <prefix>.hmm.txt and <prefix>.hmm_fract.txt. (required)"

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
            --from-bcf \
            -i ~{input_bcf} \
            -o ~{prefix} \
            --bcf-read-mode ~{bcf_read_mode} \
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
