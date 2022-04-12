version 1.0

import "Structs.wdl"

# Given BAM, call SVs using Sniffles
task Sniffles {
    input {
        File bam
        File bai
        File tandem_repeats_bed

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(1.5 * size(bam, "GB"))

    command <<<
        set -x
        num_cores=$(grep -c '^processor' /proc/cpuinfo | awk '{ print $1 - 1 }')

        sniffles -t ${num_cores} \
                 --input ~{bam} \
                 --vcf ~{prefix}.sniffles.vcf \
                 --tandem-repeats ~{tandem_repeats_bed}
    >>>

    output {
        File vcf = "~{prefix}.sniffles.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "quay.io/ymostovoy/lr-sniffles:latest"
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
