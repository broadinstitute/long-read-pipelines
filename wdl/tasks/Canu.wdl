version 1.0

import "Structs.wdl"

workflow CorrectTrimAssemble {
    input {
        String file_prefix
        String genome_size
        Array[File] reads_fastq
        Array[Float] corrected_error_rates
    }

    call Correct {
        input:
            file_prefix = file_prefix,
            genome_size = genome_size,
            reads_fastq = reads_fastq
    }

    call Trim {
        input:
            file_prefix = file_prefix,
            genome_size = genome_size,
            corrected_reads_fasta_gz = Correct.corrected_reads_fasta_gz
    }

    scatter (corrected_error_rate in corrected_error_rates) {
        call Assemble {
            input:
                genome_size = genome_size,
                file_prefix = file_prefix,
                trimmed_reads_fasta_gz = Trim.trimmed_reads_fasta_gz,
                corrected_error_rate = corrected_error_rate
        }
    }

    output {
        Array[File] assemblies = Assemble.canu_assemble_output_tgz
    }
}

task Correct {
    input {
        String file_prefix
        String genome_size
        Array[File] reads_fastq

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 30 * ceil(size(reads_fastq, "GB"))

    command <<<
        set -euxo pipefail

        /canu-2.0/Linux-amd64/bin/canu -correct \
            -p ~{file_prefix} -d canu_correct_output \
            genomeSize=~{genome_size} \
            -nanopore \
            ~{sep=' ' reads_fastq} \
            || cat /cromwell_root/monitoring.log
    >>>

    output {
        File corrected_reads_fasta_gz = "canu_correct_output/~{file_prefix}.correctedReads.fasta.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
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

task Trim {
    input {
        String file_prefix
        String genome_size
        File corrected_reads_fasta_gz

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 50 * ceil(size(corrected_reads_fasta_gz, "GB"))

    command <<<
       set -euxo pipefail

       /canu-2.0/Linux-amd64/bin/canu -trim \
             -p ~{file_prefix} -d canu_trim_output \
            genomeSize=~{genome_size} \
            -nanopore-corrected \
            ~{corrected_reads_fasta_gz} \
            || cat /cromwell_root/monitoring.log
    >>>

    output {
        File trimmed_reads_fasta_gz = "canu_trim_output/~{file_prefix}.trimmedReads.fasta.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
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

task Assemble {
    input {
        String file_prefix
        String genome_size
        File trimmed_reads_fasta_gz
        Float corrected_error_rate

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 50 * ceil(size(trimmed_reads_fasta_gz, "GB"))

    command <<<
        set -euxo pipefail

        /canu-2.0/Linux-amd64/bin/canu -assemble \
            -p ~{file_prefix} -d canu_assemble_erate_~{corrected_error_rate}_output \
            genomeSize=~{genome_size} \
            corErrorRate=~{corrected_error_rate} \
            -nanopore-corrected \
            ~{trimmed_reads_fasta_gz} \
            || cat /cromwell_root/monitoring.log

        tar czf canu_assemble_erate_~{corrected_error_rate}_output.tgz canu_assemble_erate_~{corrected_error_rate}_output/
    >>>

    output {
        File canu_assemble_output_tgz = "canu_assemble_erate_~{corrected_error_rate}_output.tgz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
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