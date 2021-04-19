version 1.0

import "Structs.wdl"

# Given BAM, call SVs using Sniffles
task Sniffles {
    input {
        File bam
        File bai

        Int min_read_support = 2
        Int min_read_length = 1000
        Int min_mq = 20

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:              "input BAM from which to call SVs"
        bai:              "index accompanying the BAM"

        min_read_support: "[default-valued] minimum reads required to make a call"
        min_read_length:  "[default-valued] filter out reads below minimum read length"
        min_mq:           "[default-valued] minimum mapping quality to accept"

        prefix:           "prefix for output"
    }

    Int cpus = 8
    Int disk_size = 4*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        sniffles -t ~{cpus} \
                 -m ~{bam} \
                 -v ~{prefix}.sniffles.vcf \
                 -s ~{min_read_support} \
                 -r ~{min_read_length} \
                 -q ~{min_mq} \
                 --num_reads_report -1 \
                 --genotype
    >>>

    output {
        File vcf = "~{prefix}.sniffles.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             15,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sv:0.1.6"
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
