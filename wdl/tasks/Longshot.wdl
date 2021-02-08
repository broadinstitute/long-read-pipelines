version 1.0

import "Structs.wdl"

# performs Longshot algo on one particular chromosome
task Longshot {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        String chr

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 4
    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB") + size(ref_fasta, "GB") + size(ref_fasta_fai, "GB"))
    String prefix = basename(bam, ".bam")

    command <<<
        set -euxo pipefail

        touch ~{prefix}.longshot.~{chr}.vcf
        longshot -F -r ~{chr} --bam ~{bam} --ref ~{ref_fasta} --out ~{prefix}.longshot.~{chr}.vcf
    >>>

    output {
        File vcf = "~{prefix}.longshot.~{chr}.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longshot:0.1.2"
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
