version 1.0

import "../../structs/Structs.wdl"

task AnnotateVCF {

    meta {
        description: "Annotate variants with contextual information"
    }

    parameter_meta {
        vcf_gz: "VCF to annotate"
        vcf_tbi: "VCF index"
        ref_fasta: "Reference sequence FASTA"
    }

    input {
        File vcf_gz
        File vcf_tbi
        File ref_fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size([vcf_gz, vcf_tbi, ref_fasta], "GB"))
    String prefix = basename(vcf_gz, ".vcf.gz") + ".annotated"

    command <<<
        set -euxo pipefail

        truvari anno gcpct -r ~{ref_fasta} ~{vcf_gz} > tmp.gcpct.vcf
        truvari anno numneigh -r 1000 -s 50 tmp.gcpct.vcf > tmp.numneigh.vcf

        cat tmp.numneigh.vcf | bgzip > ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf_gz = "~{prefix}.annotated.vcf.gz"
        File annotated_tbi = "~{prefix}.annotated.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             24,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        #docker:             "us.gcr.io/broad-dsp-lrma/lr-sv:0.1.8"
        docker:             "fcunial/truvari_intrasample:latest"
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
