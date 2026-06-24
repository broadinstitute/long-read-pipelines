version 1.0

import "../../structs/Structs.wdl"

task FilterVCFsForFlare {

    meta {
        description: "Normalize, filter, intersect, and thin study and reference VCFs for FLARE"
    }

    parameter_meta {
        gt_bcf: "Study BCF or VCF for a single chromosome"
        ref_vcf: "Reference VCF for the same chromosome"
        ref_fasta: "GRCh38 reference FASTA"
        chromosome: "Chromosome name (e.g. chr20 or chrX)"
        prefix: "Output prefix for filtered VCFs"
        is_chr_x: "Whether this chromosome is chrX"
        sample_sex_map: "Optional sample sex map for chrX handling"
        maf: "Minimum minor allele frequency in study cohort"
        thin_bp: "Minimum spacing between retained SNPs"
    }

    input {
        File gt_bcf
        File ref_vcf
        File ref_fasta
        String chromosome
        String prefix

        Boolean is_chr_x = false
        File? sample_sex_map

        Float maf = 0.01
        Int thin_bp = 20000

        RuntimeAttr? runtime_attr_override
    }

    Float gt_size_gb = size(gt_bcf, "GB")
    Float ref_size_gb = size(ref_vcf, "GB")
    Int disk_size = 4 * ceil(gt_size_gb + ref_size_gb) + 50

    String sex_map_arg = if defined(sample_sex_map) then "--sample-sex-map ${sample_sex_map}" else ""
    String chr_x_flag = if is_chr_x then "--is-chr-x" else ""

    command <<<
        set -euxo pipefail

        python3 /opt/flare/scripts/prep_vcf_for_flare.py \
            --gt ~{gt_bcf} \
            --ref ~{ref_vcf} \
            --ref-fasta ~{ref_fasta} \
            --chromosome ~{chromosome} \
            --out-prefix ~{prefix} \
            --maf ~{maf} \
            --thin-bp ~{thin_bp} \
            ~{chr_x_flag} \
            ~{sex_map_arg}
    >>>

    output {
        File gt_vcf = "~{prefix}.gt.flare.vcf.gz"
        File gt_vcf_tbi = "~{prefix}.gt.flare.vcf.gz.tbi"
        File ref_vcf_out = "~{prefix}.ref.flare.vcf.gz"
        File ref_vcf_tbi = "~{prefix}.ref.flare.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-flare:0.6.0"
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

task RunFlare {

    meta {
        description: "Run FLARE local ancestry inference for a single chromosome"
    }

    parameter_meta {
        gt_vcf: "Filtered study VCF"
        ref_vcf: "Filtered reference VCF"
        ref_panel: "FLARE ref-panel file (sample population map)"
        genetic_map: "PLINK genetic map for this chromosome"
        prefix: "Output prefix"
        em: "Whether to estimate model parameters (true for chr1, false otherwise)"
        flare_model: "Optional model file from chr1 when em=false"
        nthreads: "Number of FLARE threads"
        mem_gb: "Java heap size in GB"
    }

    input {
        File gt_vcf
        File ref_vcf
        File ref_panel
        File genetic_map
        String prefix

        Boolean em
        File? flare_model

        Int nthreads = 16
        Int mem_gb = 128

        RuntimeAttr? runtime_attr_override
    }

    Float input_size_gb = size(gt_vcf, "GB") + size(ref_vcf, "GB")
    Int disk_size = 4 * ceil(input_size_gb) + 50

    String model_arg = if defined(flare_model) then "model=~{flare_model}" else ""

    command <<<
        set -euxo pipefail

        java -Xmx~{mem_gb}g -jar /opt/flare/flare.jar \
            ref=~{ref_vcf} \
            gt=~{gt_vcf} \
            map=~{genetic_map} \
            ref-panel=~{ref_panel} \
            out=~{prefix} \
            nthreads=~{nthreads} \
            em=~{em} \
            ~{model_arg}
    >>>

    output {
        File global_anc = "~{prefix}.global.anc.gz"
        File model = "~{prefix}.model"
        File anc_vcf = "~{prefix}.anc.vcf.gz"
        File log = "~{prefix}.log"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          nthreads,
        mem_gb:             mem_gb,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-flare:0.6.0"
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
