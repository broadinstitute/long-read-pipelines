version 1.0

import "../../structs/Structs.wdl"


workflow DeepVariant {

    meta {
        description: "Workflow for getting VCF and gVCF from DeepVariant. Note VCF is un-phased."
    }

    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        Int pepper_threads
        Int pepper_memory

        Int dv_threads
        Int dv_memory

        String zones = "us-central1-b us-central1-c"
    }

    parameter_meta {
        zones: "select which zone (GCP) to run this task"
    }

    call DV as deep_variant {
        input:
            bam = bam,
            bai = bai,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            threads = dv_threads,
            memory = dv_memory,
            zones = zones
    }

    output {
        File VCF        = deep_variant.VCF
        File VCF_tbi    = deep_variant.VCF_tbi

        File gVCF       = deep_variant.gVCF
        File gVCF_tbi   = deep_variant.gVCF_tbi
    }
}

task DV {

    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        Int threads
        Int memory
        String zones

        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, ".bam") + ".deepvariant"
    String output_root = "/cromwell_root/dv_output"

    Int bam_sz = ceil(size(bam, "GB"))
    Boolean is_big_bam = bam_sz > 100
    Int inflation_factor = if (is_big_bam) then 10 else 5
    Int minimal_disk = 1000
	Int disk_size = if inflation_factor * bam_sz > minimal_disk then inflation_factor * bam_sz else minimal_disk

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        mkdir -p "~{output_root}"

        /opt/deepvariant/bin/run_deepvariant \
            --model_type=WGS \
            --ref=~{ref_fasta} \
            --reads=~{bam} \
            --output_vcf="~{output_root}/~{prefix}.vcf.gz" \
            --output_gvcf="~{output_root}/~{prefix}.g.vcf.gz" \
            --num_shards="${num_core}"

        find "~{output_root}/" -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g' \
            > "~{output_root}/dir_structure.txt"
    >>>

    output {
        File output_dir_structure = "~{output_root}/dir_structure.txt"

        File VCF        = "~{output_root}/~{prefix}.vcf.gz"
        File VCF_tbi    = "~{output_root}/~{prefix}.vcf.gz.tbi"

        File gVCF       = "~{output_root}/~{prefix}.g.vcf.gz"
        File gVCF_tbi   = "~{output_root}/~{prefix}.g.vcf.gz.tbi"

        File visual_report_html = "~{output_root}/~{prefix}.visual_report.html"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          threads,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "google/deepvariant:1.4.0"
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


task DeepTrio {

    input {
        File parent1_bam
        File parent1_bai
        String parent1_sample_id

        File parent2_bam
        File parent2_bai
        String parent2_sample_id

        File proband_bam
        File proband_bai
        String proband_sample_id

        File ref_fasta
        File ref_fasta_fai
        File ref_fasta_dict

        String model_type = "WGS"

        RuntimeAttr? runtime_attr_override
    }

    String output_root = "/cromwell_root/dv_output"

    # Calculate a bunch of space:
	Int disk_size = 10 + (5 * ceil(ceil(size(parent1_bam, "GB")) + ceil(size(parent2_bam, "GB")) + ceil(size(proband_bam, "GB")) + ceil(size(ref_fasta, "GB"))))

    command <<<
        set -euxo pipefail

        num_core=$(awk '/^processor/{print $3}' /proc/cpuinfo | wc -l)

        mkdir -p "~{output_root}/logging"

        /opt/deepvariant/bin/deeptrio/run_deeptrio \
            --model_type=~{model_type} \
            --num_shards "${num_core}"  \
            --logging_dir ~{output_root}/logging \
            --intermediate_results_dir ~{output_root}/intermediate_results_dir \
            \
            --ref=~{ref_fasta} \
            \
            --reads_child=~{proband_bam} \
            --reads_parent1=~{parent1_bam} \
            --reads_parent2=~{parent2_bam} \
            \
            --sample_name_child "~{proband_sample_id}" \
            --output_vcf_child ~{output_root}/~{proband_sample_id}.output.vcf.gz \
            --output_gvcf_child ~{output_root}/~{proband_sample_id}.g.vcf.gz \
            \
            --sample_name_parent1 "~{parent1_sample_id}" \
            --output_vcf_parent1 ~{output_root}/~{parent1_sample_id}.output.vcf.gz \
            --output_gvcf_parent1 ~{output_root}/~{parent1_sample_id}.g.vcf.gz \
            \
            --sample_name_parent2 "~{parent2_sample_id}" \
            --output_vcf_parent2 ~{output_root}/~{parent2_sample_id}.output.vcf.gz \
            --output_gvcf_parent2 ~{output_root}/~{parent2_sample_id}.g.vcf.gz \
            \
            --dry_run=false \

        cd ~{output_root}
        tar -zcf logging.tar.gz logging
    >>>

    output {
        File parent1_vcf        = "~{output_root}/~{parent1_sample_id}.output.vcf.gz"
        File parent1_vcf_index  = "~{output_root}/~{parent1_sample_id}.output.vcf.gz.tbi"
        File parent1_gvcf       = "~{output_root}/~{parent1_sample_id}.output.g.vcf.gz"
        File parent1_gvcf_index = "~{output_root}/~{parent1_sample_id}.output.g.vcf.gz.tbi"

        File parent2_vcf        = "~{output_root}/~{parent2_sample_id}.output.vcf.gz"
        File parent2_vcf_index  = "~{output_root}/~{parent2_sample_id}.output.vcf.gz.tbi"
        File parent2_gvcf       = "~{output_root}/~{parent2_sample_id}.output.g.vcf.gz"
        File parent2_gvcf_index = "~{output_root}/~{parent2_sample_id}.output.g.vcf.gz.tbi"

        File proband_vcf        = "~{output_root}/~{proband_sample_id}.output.vcf.gz"
        File proband_vcf_index  = "~{output_root}/~{proband_sample_id}.output.vcf.gz.tbi"
        File proband_gvcf       = "~{output_root}/~{proband_sample_id}.output.g.vcf.gz"
        File proband_gvcf_index = "~{output_root}/~{proband_sample_id}.output.g.vcf.gz.tbi"

        File logs = "~{output_root}/logging.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "google/deepvariant:deeptrio-1.8.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
