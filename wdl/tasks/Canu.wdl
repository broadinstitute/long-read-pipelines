version 1.0

import "Structs.wdl"

task CorrectTrimAssemble {
    input {
        String file_prefix
        String genome_size
        Array[File] reads_fastq
        Float corrected_error_rate

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 50 * ceil(size(reads_fastq, "GB"))

    command <<<
        set -euxo pipefail

        /canu-2.0/Linux-amd64/bin/canu -correct \
            -p ~{file_prefix} -d ~{file_prefix}_correct \
            genomeSize=~{genome_size} \
            -nanopore \
            ~{sep=' ' reads_fastq} \
            || cat /cromwell_root/monitoring.log

        /canu-2.0/Linux-amd64/bin/canu -trim \
            -p ~{file_prefix} -d ~{file_prefix}_trim \
            genomeSize=~{genome_size} \
            -nanopore-corrected \
            ~{file_prefix}_correct/~{file_prefix}.correctedReads.fasta.gz \
            || cat /cromwell_root/monitoring.log

        /canu-2.0/Linux-amd64/bin/canu -assemble \
            -p ~{file_prefix} -d ~{file_prefix}_assemble \
            genomeSize=~{genome_size} \
            corErrorRate=~{corrected_error_rate} \
            -nanopore-corrected \
            ~{file_prefix}_trim/~{file_prefix}.trimmedReads.fasta.gz \
            || cat /cromwell_root/monitoring.log
    >>>

    output {
        File assembly = "~{file_prefix}_assemble/~{file_prefix}.contigs.fasta"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-canu:0.1.0"
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


