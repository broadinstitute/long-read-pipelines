version 1.0

import "Structs.wdl"

task Assemble {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 * ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools fasta ~{bam} > reads.fasta
        shasta-Linux-0.2.0 --input reads.fasta

        find . -type f -exec ls -lah '{}' ';'
    >>>

    output {
        String out = read_lines(stdout())
        File assembled_fasta = "./ShastaRun/Assembly.fasta"
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:      8,
        mem_gb:         64,
        disk_gb:        disk_size,
        boot_disk_gb:   10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-asm:0.01.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:        select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:     select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:     select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:        select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:         select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:             select_first([runtime_attr.docker, default_attr.docker])
    }
}