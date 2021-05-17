version 1.0

import "Structs.wdl"
import "Utils.wdl" as Utils

workflow C3POa {
    input {
        File manifest_chunk
        File ref_fasta
        File splint_fasta
    }

    call Cat as CatRawReads { input: files = read_lines(manifest_chunk), out = "chunk.fastq" }

    call Processing { input: fastq = CatRawReads.merged, splint_fasta = splint_fasta }

    call Cat as CatSubreads1 { input: files = Processing.subreads1, out = "subreads1.fastq" }
    call Cat as CatSubreads2 { input: files = Processing.subreads2, out = "subreads2.fastq" }
    call Cat as CatSubreads3 { input: files = Processing.subreads3, out = "subreads3.fastq" }
    call Cat as CatSubreads4 { input: files = Processing.subreads4, out = "subreads4.fastq" }

    call Cat as CatConsensus1 { input: files = Processing.consensus1, out = "consensus1.fasta" }
    call Cat as CatConsensus2 { input: files = Processing.consensus2, out = "consensus2.fasta" }
    call Cat as CatConsensus3 { input: files = Processing.consensus3, out = "consensus3.fasta" }
    call Cat as CatConsensus4 { input: files = Processing.consensus4, out = "consensus4.fasta" }

    output {
        Array[File] subreads1 = CatSubreads1.merged
        Array[File] subreads2 = CatSubreads2.merged
        Array[File] subreads3 = CatSubreads3.merged
        Array[File] subreads4 = CatSubreads4.merged

        Array[File] consensus1 = CatConsensus1.merged
        Array[File] consensus2 = CatConsensus2.merged
        Array[File] consensus3 = CatConsensus3.merged
        Array[File] consensus4 = CatConsensus4.merged
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
            -l 1000 -d 500 -n 32 -g 1000 \
            -o out

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
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
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
        boot_disk_gb:       10,
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
        boot_disk_gb:       10,
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
