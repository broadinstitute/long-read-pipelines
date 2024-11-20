version 1.0

import "../../structs/Structs.wdl"
import "../../tasks/Utility/Utils.wdl" as Utils

workflow C3POa {
    input {
        File manifest_chunk
        File ref_fasta
        File splint_fasta
    }

    call Cat as CatRawReads { input: files = read_lines(manifest_chunk), out = "chunk.fastq" }

    call Processing { input: fastq = CatRawReads.merged, splint_fasta = splint_fasta }

    output {
        File subreads1 = Processing.subreads1
        File subreads2 = Processing.subreads2
        File subreads3 = Processing.subreads3
        File subreads4 = Processing.subreads4

        File consensus1 = Processing.consensus1
        File consensus2 = Processing.consensus2
        File consensus3 = Processing.consensus3
        File consensus4 = Processing.consensus4

        Int no_splint_reads  = Processing.no_splint_reads
        Int under_len_cutoff = Processing.under_len_cutoff
        Int total_reads      = Processing.total_reads
    }
}

task Processing {
    input {
        File fastq
        File splint_fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(fastq, "GB"))

    command <<<
        set -euxo pipefail

        mkdir out
        python3 /C3POa/C3POa.py \
            -r ~{fastq} \
            -s ~{splint_fasta} \
            -c /c3poa.config.txt \
            -l 100 -d 500 -n 32 -g 1000 \
            -o out

        grep 'No splint reads' out/c3poa.log | awk '{ print $4 }' > no_splint_reads.txt
        grep 'Under len cutoff' out/c3poa.log | awk '{ print $4 }' > under_len_cutoff.txt
        grep 'Total reads' out/c3poa.log | awk '{ print $3 }' > total_reads.txt

        tree -h
    >>>

    output {
        File consensus1 = "out/10x_Splint_1/R2C2_Consensus.fasta"
        File consensus2 = "out/10x_Splint_2/R2C2_Consensus.fasta"
        File consensus3 = "out/10x_Splint_3/R2C2_Consensus.fasta"
        File consensus4 = "out/10x_Splint_4/R2C2_Consensus.fasta"

        File subreads1 = "out/10x_Splint_1/R2C2_Subreads.fastq"
        File subreads2 = "out/10x_Splint_2/R2C2_Subreads.fastq"
        File subreads3 = "out/10x_Splint_3/R2C2_Subreads.fastq"
        File subreads4 = "out/10x_Splint_4/R2C2_Subreads.fastq"

        File c3poa_log = "out/c3poa.log"
        Int no_splint_reads = read_int("no_splint_reads.txt")
        Int under_len_cutoff = read_int("under_len_cutoff.txt")
        Int total_reads = read_int("total_reads.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-c3poa:2.2.2"
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

        python3 /C3POa/C3POa_postprocessing.py -i ~{consensus} -c /c3poa.config.txt -a /C3POa/adapter.fasta -o ./
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
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-c3poa:0.1.9"
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

    Int disk_size = 1 + 3*ceil(size(files, "GB"))

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
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-c3poa:0.1.9"
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
