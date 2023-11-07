version 1.0

import "../../structs/Structs.wdl"

task VariantFiltration {
    meta {
        desciption: "A simple task to apply a filter annotation (to the FILTER column) to an VCF"
    }
    parameter_meta {
        vcf: "vcf to be annotated"
        out_appendix: "appendix to add to the input VCF's basename"
        filter_exp: "argument to be passed to `--filter-expression` of VariantFiltration in GATK"
        filter_name: "argument to be passed to `--filter-name` of VariantFiltration in GATK"
        gatk_docker_tag: "Tag of the GATK GCR image, hosted at us.gcr.io/broad-gatk/gatk. Example: 4.4.0.0"
        annotated_vcf: "the annotated VCF"
    }
    input {
        File vcf
        File tbi
        File ref_fasta
        File ref_fai
        File ref_dict
        String out_appendix
        String filter_exp
        String filter_name
        String gatk_docker_tag

        RuntimeAttr? runtime_attr_override
    }
    output {
        File  annotated_vcf = out
        File? annotated_tbi = "~{out}.tbi"
    }

    String base = basename(vcf, ".vcf.gz")
    String out = "~{base}.~{out_appendix}.vcf.gz"

    command <<<
        set -euxo pipefail

        gatk \
        VariantFiltration \
            -R ~{ref_fasta} \
            --output "~{out}" \
            --variant "~{vcf}" \
            --filter-expression "~{filter_exp}" \
            --filter-name "~{filter_name}"
    >>>

    #########################
    Int disk_size = 10 + 4*ceil(size(vcf, "GiB"))
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-gatk/gatk:~{gatk_docker_tag}"
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