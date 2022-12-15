version 1.0

import "Structs.wdl"

task BamToFq {
    input {
        File bam
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools sort -n ~{bam} | samtools bam2fq \
            -n \
            -s /dev/null \
            -1 ~{prefix}.end1.fq.gz \
            -2 ~{prefix}.end2.fq.gz \
            -0 ~{prefix}.unpaired.fq.gz
    >>>

    output {
        File fq_end1 = "~{prefix}.end1.fq.gz"
        File fq_end2 = "~{prefix}.end1.fq.gz"
        File fq_unpaired = "~{prefix}.unpaired.fq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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
