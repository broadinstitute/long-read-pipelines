version 1.0

import "Structs.wdl"

workflow Process {
    input {
        File input_fasta
        Int parallel_instances
    }

    call SplitReads {
        input:
            combined_read_fasta = input_fasta,
            parallel_instances = parallel_instances
    }

    scatter (split_read_fasta in SplitReads.split_reads_fasta) {
        call RemovePalindromes {
            input:
                read_fasta = split_read_fasta
        }
    }

    call MergeFasta {
        input:
            processed_fastas = RemovePalindromes.processed_fasta
    }

    output {
        File processed_fasta = MergeFasta.processed_fasta
    }
}

task SplitReads {
    input {
        File combined_read_fasta
        Int parallel_instances
    }

    Int disk_size = 3 * ceil(size(combined_read_fasta, "GB"))

    command <<<
        set -euxo pipefail

        seqkit split2 ~{combined_read_fasta} -p ~{parallel_instances} -O splits
    >>>

    output {
        Array[File] split_reads_fasta = glob("splits/*")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0
    }

    runtime {
        cpu:                    default_attr.cpu_cores
        memory:                 default_attr.mem_gb + " GiB"
        disks: "local-disk " +  default_attr.disk_gb + " HDD"
        bootDiskSizeGb:         default_attr.boot_disk_gb
        preemptible:            default_attr.preemptible_tries
        maxRetries:             default_attr.max_retries
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
        mem_gb:             3,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  3,
        max_retries:        0
    }

    runtime {
        cpu:                    default_attr.cpu_cores
        memory:                 default_attr.mem_gb + " GiB"
        disks: "local-disk " +  default_attr.disk_gb + " HDD"
        bootDiskSizeGb:         default_attr.boot_disk_gb
        preemptible:            default_attr.preemptible_tries
        maxRetries:             default_attr.max_retries
        gpuType:                "nvidia-tesla-t4"
        gpuCount:               1
        nvidiaDriverVersion:    "418.87.00"
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

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0
    }
    
    runtime {
        cpu:                    default_attr.cpu_cores
        memory:                 default_attr.mem_gb + " GiB"
        disks: "local-disk " +  default_attr.disk_gb + " HDD"
        bootDiskSizeGb:         default_attr.boot_disk_gb
        preemptible:            default_attr.preemptible_tries
        maxRetries:             default_attr.max_retries
        docker:                 "quay.io/broad-long-read-pipelines/lr-pacasus:0.3.0"
    }
}