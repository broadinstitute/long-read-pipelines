version 1.0

##########################################################################################
# A workflow that runs the Canu 3-step assembly (correct, trim, assemble).
# - Tested on a small genome (malaria ~23mb), larger genomes will likely require some changes
# including tweaks to the default resource allocation.
# - Currently assumes nanopore reads
##########################################################################################

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.36/wdl/tasks/Structs.wdl"

workflow CorrectTrimAssemble {
    input {
        String output_file_prefix
        String genome_size
        File reads
        Float correct_corrected_error_rate
        Float trim_corrected_error_rate
        Float assemble_corrected_error_rate
    }

    call Correct {
        input:
            output_file_prefix = output_file_prefix,
            genome_size = genome_size,
            reads = reads,
            corrected_error_rate = correct_corrected_error_rate
    }

    call Trim {
        input:
            output_file_prefix = output_file_prefix,
            genome_size = genome_size,
            corrected_reads_fasta_gz = Correct.corrected_reads_fasta_gz,
            corrected_error_rate = trim_corrected_error_rate
    }

    call Assemble {
        input:
            genome_size = genome_size,
            output_file_prefix = output_file_prefix,
            trimmed_reads_fasta_gz = Trim.trimmed_reads_fasta_gz,
            corrected_error_rate = assemble_corrected_error_rate
    }

    output {
        File canu_contigs_fasta = Assemble.canu_contigs_fasta
    }
}

# performs canu correct on raw reads, currently assumes ONT reads
task Correct {
    input {
        String output_file_prefix
        String genome_size
        File reads
        Float corrected_error_rate

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:                  "reads to be canu-corrected"
        genome_size:            "estimate on genome size (parameter to canu's \'genomeSize\')"
        output_file_prefix:     "prefix to output files"
        corrected_error_rate:   "parameter to canu's \'correctedErrorRate\'"
    }

    Int disk_size = 100 * ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        canu -correct \
            -p ~{output_file_prefix} -d canu_correct_output \
            genomeSize=~{genome_size} \
            corMaxEvidenceErate=0.15 \
            correctedErrorRate=~{corrected_error_rate} \
            -nanopore \
            ~{reads} \
            || cat /cromwell_root/monitoring.log
    >>>

    output {
        File corrected_reads_fasta_gz = "canu_correct_output/~{output_file_prefix}.correctedReads.fasta.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-canu:0.1.0"
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

# performs canu trim on corrected reads
task Trim {
    input {
        String output_file_prefix
        String genome_size
        File corrected_reads_fasta_gz
        Float corrected_error_rate

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        corrected_reads_fasta_gz:   "reads that have been canu-corrected"
        genome_size:                "estimate on genome size (parameter to canu's \'genomeSize\')"
        output_file_prefix:         "prefix to output files"
        corrected_error_rate:       "parameter to canu's \'correctedErrorRate\'"
    }

    Int disk_size = 50 * ceil(size(corrected_reads_fasta_gz, "GB"))

    command <<<
       set -euxo pipefail

       canu -trim \
             -p ~{output_file_prefix} -d canu_trim_output \
            genomeSize=~{genome_size} \
            correctedErrorRate=~{corrected_error_rate} \
            -nanopore-corrected \
            ~{corrected_reads_fasta_gz} \
            || cat /cromwell_root/monitoring.log
    >>>

    output {
        File trimmed_reads_fasta_gz = "canu_trim_output/~{output_file_prefix}.trimmedReads.fasta.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-canu:0.1.0"
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

# performs assembly on corrected, then trimmmed reads
task Assemble {
    input {
        String output_file_prefix
        String genome_size
        File trimmed_reads_fasta_gz
        Float corrected_error_rate

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        trimmed_reads_fasta_gz:     "reads that have been canu-corrected-trimmed"
        genome_size:                "estimate on genome size (parameter to canu's \'genomeSize\')"
        output_file_prefix:         "prefix to output files"
        corrected_error_rate:       "parameter to canu's \'correctedErrorRate\'"
    }

    Int disk_size = 50 * ceil(size(trimmed_reads_fasta_gz, "GB"))

    command <<<
        set -euxo pipefail

        canu -assemble \
            -p ~{output_file_prefix} -d canu_assemble_output \
            genomeSize=~{genome_size} \
            correctedErrorRate=~{corrected_error_rate} \
            -nanopore-corrected \
            ~{trimmed_reads_fasta_gz} \
            || cat /cromwell_root/monitoring.log
    >>>

    output {
        File canu_contigs_fasta = "canu_assemble_output/~{output_file_prefix}.contigs.fasta"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-canu:0.1.0"
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
