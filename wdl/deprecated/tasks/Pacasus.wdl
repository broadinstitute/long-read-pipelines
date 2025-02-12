version 1.0

import "../../tasks/Utility/Utils.wdl" as Utils
import "../../structs/Structs.wdl"

workflow Process {
    input {
        File reads
        Int chunk_size_mb
    }

    call Utils.ConvertReads {
        input:
            reads = reads,
            output_format = "fasta"
    }

    call SplitReadsByLength {
        input:
            reads_fasta = ConvertReads.converted_reads,
            length = 55000
    }

    call SplitReads {
        input:
            reads_fasta = SplitReadsByLength.shorter_reads,
            chunk_size_mb = chunk_size_mb
    }

    scatter (split_read_fasta in SplitReads.split_reads_fasta) {
        call RemovePalindromes {
            input:
                read_fasta = split_read_fasta
        }
    }

    call MergeFasta {
        input:
            processed_fastas = flatten([RemovePalindromes.processed_fasta, [SplitReadsByLength.longer_reads]])
    }

    output {
        File processed_fasta = MergeFasta.processed_fasta
    }
}

task SplitReadsByLength {
    input {
        File reads_fasta
        Int length
    }

    Int disk_size = 3 * ceil(size(reads_fasta, "GB"))

    command <<<
        set -euxo pipefail

        seqkit seq ~{reads_fasta} --max-len $(expr ~{length} + 1) > shorter_reads.fasta
        seqkit seq ~{reads_fasta} --min-len ~{length} > longer_reads.fasta
    >>>

    output {
        File shorter_reads = "shorter_reads.fasta"
        File longer_reads = "longer_reads.fasta"
    }

    runtime {
        cpu:                    4
        memory:                 "8 GiB"
        disks:                  "local-disk " +  disk_size + " HDD"
        bootDiskSizeGb:         25
        preemptible:            2
        maxRetries:             0
        docker:                 "quay.io/broad-long-read-pipelines/lr-pacasus:0.3.0"
    }
}

task SplitReads {
    input {
        File reads_fasta
        Int chunk_size_mb
    }

    Int mem_size = 2 * ceil(size(reads_fasta, "GB"))
    Int disk_size = 4 * ceil(size(reads_fasta, "GB"))

    command <<<
        set -euxo pipefail

        reads_size_mb=$(ls -s --b=M ~{reads_fasta} | cut -d'M' -f1)
        n_chunks=$(expr $reads_size_mb / ~{chunk_size_mb})
        cat ~{reads_fasta} | seqkit split -p $n_chunks -O splits
    >>>

    output {
        Array[File] split_reads_fasta = glob("splits/*")
    }

    runtime {
        cpu:                    4
        memory:                 mem_size + " GiB"
        disks:                  "local-disk " +  disk_size + " HDD"
        bootDiskSizeGb:         25
        preemptible:            0
        maxRetries:             0
        docker:                 "quay.io/broad-long-read-pipelines/lr-pacasus:0.3.0"
    }
}

task RemovePalindromes {
    input {
        File read_fasta
    }

    Int disk_size = 3 * ceil(size(read_fasta, "GB"))
    String read_basename = basename(read_fasta)

    command <<<
        set -euxo pipefail

        python /pacasus/pacasus.py --device_type=GPU --platform_name=NVIDIA --framework=cuda ~{read_fasta} -o ~{read_basename}.processed.fasta
    >>>

    output {
        File processed_fasta = "~{read_basename}.processed.fasta"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             6,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0
    }

    runtime {
        cpu:                    default_attr.cpu_cores
        memory:                 select_first([default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         default_attr.boot_disk_gb
        preemptible:            default_attr.preemptible_tries
        maxRetries:             default_attr.max_retries
        gpuType:                "nvidia-tesla-t4"
        gpuCount:               1
        nvidiaDriverVersion:    "418.152.00"
        zones:                  ["us-east1-c"]
        docker:                 "quay.io/broad-long-read-pipelines/lr-pacasus:0.3.0"
    }
}

task MergeFasta {
    input {
        Array[File] processed_fastas
    }

    Int disk_size = 3 * ceil(size(processed_fastas, "GB"))

    command <<<
        set -euxo pipefail

        cat ~{sep=" " processed_fastas} > pacasus_processed.fasta
    >>>

    output {
       File processed_fasta = "pacasus_processed.fasta"
    }

    runtime {
        cpu:                    4
        memory:                 "8 GiB"
        disks:                  "local-disk " +  disk_size + " HDD"
        bootDiskSizeGb:         25
        preemptible:            2
        maxRetries:             0
        docker:                 "quay.io/broad-long-read-pipelines/lr-pacasus:0.3.0"
    }
}

