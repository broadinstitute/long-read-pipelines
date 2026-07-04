version 1.0

import "../../structs/Structs.wdl"


task FilterVCFsForFlare {

    input {
        File joint_vcf
        File ref_vcf
        File ref_fasta
        File ref_fasta_fai
        String chromosome
        String prefix

        Boolean is_chr_x = false
        File? sample_sex_map

        Float maf = 0.01
        Int thin_bp = 20000
    }

    Float gt_size_gb = size(joint_vcf, "GB")
    Float ref_size_gb = size(ref_vcf, "GB")
    Int disk_size = 4 * ceil(gt_size_gb + ref_size_gb) + 50

    String sex_map_arg = if defined(sample_sex_map) then "--sample-sex-map ${sample_sex_map}" else ""
    String chr_x_flag = if is_chr_x then "--is-chr-x" else ""

    command <<<
        set -euxo pipefail

        python3 /opt/flare/scripts/prep_vcf_for_flare.py \
            --gt ~{joint_vcf} \
            --ref ~{ref_vcf} \
            --ref-fasta ~{ref_fasta} \
            --ref-fasta-fai ~{ref_fasta_fai} \
            --chromosome ~{chromosome} \
            --out-prefix ~{prefix} \
            --maf ~{maf} \
            --thin-bp ~{thin_bp} \
            ~{chr_x_flag} \
            ~{sex_map_arg}
    >>>

    output {
        File gt_vcf = "~{prefix}.gt.flare.vcf.gz"
        File gt_vcf_csi = "~{prefix}.gt.flare.vcf.gz.csi"
        File ref_vcf_out = "~{prefix}.ref.flare.vcf.gz"
        File ref_vcf_csi = "~{prefix}.ref.flare.vcf.gz.csi"
    }

    runtime {
        cpu: 8
        memory: "64 GiB"
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: 10
        preemptible: 1
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-flare:0.6.0-v3"
    }
}

task Flare {

    meta {
        description: "Local Ancestry Inference"
    }


    input {
        File ref_vcf
        File ref_vcf_index
        File ref_panel
        File test_vcf
        File test_vcf_index
        File plink_map
        String output_prefix
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        Boolean em = true
        File? flare_model

        Int nthreads = 16
        Int mem_gb = 64

        RuntimeAttr? runtime_attr_override
    }

    String model_arg = if defined(flare_model) then "model=~{flare_model}" else ""

    command <<<
        set -euxo pipefail

        java -Xmx~{mem_gb}g -jar /LAI/flare.jar \
            ref=~{ref_vcf} \
            gt=~{test_vcf} \
            map=~{plink_map} \
            ref-panel=~{ref_panel} \
            out=~{output_prefix} \
            nthreads=~{nthreads} \
            em=~{em} \
            ~{model_arg}

    >>>

    output {
        File global_anc = "~{output_prefix}.global.anc.gz"
        File model = "~{output_prefix}.model"
        File anc_vcf = "~{output_prefix}.anc.vcf.gz"
        File log = "~{output_prefix}.log"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          nthreads,
        mem_gb:             mem_gb,
        disk_gb:            200,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "hangsuunc/flare:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
