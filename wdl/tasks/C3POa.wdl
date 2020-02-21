version 1.0

import "Structs.wdl"
import "Utils.wdl" as Utils

workflow C3POa {
    input {
        File manifest_chunk
        File ref_fasta
    }

    call Preprocessing { input: manifest_chunk = manifest_chunk }

    scatter (fastq in Preprocessing.fastqs) {
        call Processing { input: preprocessed_fastq = fastq }
    }

    call Cat as CatSubreads { input: files = Processing.subreads, out = "subreads.fastq" }
    call Cat as CatConsensus { input: files = Processing.consensus, out = "consensus.fasta" }

    output {
        File subreads = CatSubreads.merged
        File consensus = CatConsensus.merged
    }
}

task Preprocessing {
    input {
        File manifest_chunk

        RuntimeAttr? runtime_attr_override
    }

    Array[String] fastqs = read_lines(manifest_chunk)

    Int disk_size = 1000

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        gsutil cat ~{sep=' ' fastqs} > chunk.fastq

        mkdir pre

        python3 /C3POa/C3POa_preprocessing.py -i chunk.fastq \
                                              -o pre \
                                              -q 9 \
                                              -l 1000 \
                                              -s /C3POa/splint.fasta \
                                              -c /c3poa.config.txt

        find pre -name 'R2C2_raw_reads.fastq' | awk -F"/" '{ print $0, "pre_" $2 "_" $4 }' | xargs -L 1 cp
    >>>

    output {
        Array[File] fastqs = glob("pre_*.fastq")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-c3poa:0.01.03"
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

task Processing {
    input {
        File preprocessed_fastq

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(preprocessed_fastq, "GB"))

    command <<<
        set -euxo pipefail

        python3 /C3POa/C3POa.py -z -d 500 -l 1000 -p ./ \
                                -r ~{preprocessed_fastq} \
                                -m /C3POa/NUC.4.4.mat \
                                -c /c3poa.config.txt \
                                -o consensus.fasta
    >>>

    output {
        File consensus = "consensus.fasta"
        File subreads = "subreads.fastq"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-c3poa:0.01.03"
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

task Postprocessing {
    input {
        File consensus

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(consensus, "GB"))

    command <<<
        set -euxo pipefail

        df -h .
        tree -h .

        python3 /C3POa/C3POa_postprocessing.py -i ~{consensus} -c /c3poa.config.txt -a /C3POa/adapter.fasta -o ./

        df -h .
        tree -h .
    >>>

    output {
        File consensus_full = "R2C2_full_length_consensus_reads.fasta"
        File consensus_left = "R2C2_full_length_consensus_reads_left_splint.fasta"
        File consensus_right = "R2C2_full_length_consensus_reads_right_splint.fasta"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-c3poa:0.01.03"
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

task Cat {
    input {
        Array[File] files
        String out

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3*ceil(size(files, "GB"))

    command <<<
        set -euxo pipefail

        cat ~{sep=' ' files} > ~{out}
    >>>

    output {
        File merged = "~{out}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-c3poa:0.01.03"
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
