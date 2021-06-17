version 1.0

import "Structs.wdl"

task Quantify {
    input {
        File aligned_bam
        File aligned_bai
        File gtf
        Boolean keep_retained_introns = false
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size([aligned_bam, aligned_bai, gtf], "GB"))

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        stringtie \
            -Lv -p $num_core ~{true='-i' false='' keep_retained_introns} \
            -G ~{gtf} \
            -o ~{prefix}.gtf \
            -A ~{prefix}.gene_abund.out \
            ~{aligned_bam}
    >>>

    output {
        File st_gtf = "~{prefix}.gtf"
        File st_abund = "~{prefix}.gene_abund.out"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-stringtie2:2.1.6"
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