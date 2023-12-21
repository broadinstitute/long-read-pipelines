version 1.0

import "../../structs/Structs.wdl"


task Stats {
    meta {
        desription:
        "seqkit stats command"
    }
    parameter_meta {
        fastq: "file to collect stats on"
        seq_type: "argument to the --seq-type paramter"
    }

    input {
        File fastq
        String seq_type
        RuntimeAttr? runtime_attr_override
    }

    output {
        Map[String, Float] res = read_map("2.col.map.tsv")
    }

    Int disk_size = 10 + ceil(size(fastq, "GB"))

    command <<<
        set -eux

        seqkit stats \
            -aT \
            -t ~{seq_type} \
            -o 2_line.tsv \
            ~{fastq}

        datamash transpose \
            < 2_line.tsv \
            | grep -vw "^file" \
            | grep -vw "^format" \
            | grep -vw "^type" \
        > 2.col.map.tsv
    >>>

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:             4,
        mem_gb:                4,
        disk_gb:               disk_size,
        preemptible_tries:     2,
        max_retries:           0,
        docker:                "us.gcr.io/broad-dsp-lrma/lr-seqkit:2.4.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                   select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        preemptible:           select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:            select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                select_first([runtime_attr.docker, default_attr.docker])
    }
}

task FilterByLength {
    meta {
        desciption: "Filter FASTQ by a length threshold (>=)."
    }
    parameter_meta {
        threshold: "Sequences shorter than this will be dropped."
        res: "result file"
    }
    input {
        File fq
        Int threshold
    }
    output {
        File res = "~{prefix}.length-filter-ge-~{threshold}.fq.gz"
    }

    String prefix = basename(basename(basename(fq, ".gz"), ".fastq"), ".fq")
    Int disk_space = 10 + 2 * ceil(size(fq, "GiB"))
    command <<<
        set -eux

        seqkit seq -m ~{threshold} ~{fq} \
        | gzip \
        > "~{prefix}.length-filter-ge-~{threshold}.fq.gz"
    >>>
    runtime {
        cpu:            2
        memory:         "8 GiB"
        disks:          "local-disk ~{disk_space} SSD"
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-seqkit:2.4.0"
    }
}

task FilterByLenSeqTk {
    meta {
        desciption:
        "Alternative implementation to FilterByLength, here using seqtk"
    }
    parameter_meta {
        exclude_len_threshold: "Sequeces shorter than this will be dropped from analysis."
    }

    input {
        File fastq
        Int exclude_len_threshold
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + 2*ceil(size(fastq, "GB"))

    String base = basename(basename(fastq, ".fastq.gz"), ".fq.gz")
    String out_prefx = "~{base}.RL_gt_~{exclude_len_threshold}"

    command <<<
        set -eux

        seqtk seq \
            -L ~{exclude_len_threshold} \
            ~{fastq} \
        | gzip \
        > "~{out_prefx}.fastq.gz"
    >>>

    output {
        File res = "~{out_prefx}.fastq.gz"
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:             4,
        mem_gb:                16,
        disk_gb:               disk_size,
        preemptible_tries:     0,
        max_retries:           0,
        docker:                "us.gcr.io/broad-dsp-lrma/lr-seqtk:1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                   select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        preemptible:           select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:            select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                select_first([runtime_attr.docker, default_attr.docker])
    }
}

task MergeFastqs {

    meta {
        description : "Merge fastq files."
    }

    parameter_meta {
        fastqs: {
            desciption: "Fastq files to be merged.",
            localization_optional: true
        }
        prefix: "Prefix for the output fastq file."
        disk_type: "type of disk for the VM"
    }

    input {
        Array[File] fastqs
        String prefix

        String disk_type
        RuntimeAttr? runtime_attr_override
    }

    command <<<
    set -eux

        mkdir -p input_fastqs
        time \
        gcloud storage cp \
            ~{sep=' ' fastqs} \
            input_fastqs/

        cd input_fastqs && ls
        FILE="~{fastqs[0]}"
        if [[ "$FILE" =~ \.gz$ ]]; then
            cat * > ~{prefix}.fq.gz
        else
            cat * | pigz > ~{prefix}.fq.gz
        fi
        mv ~{prefix}.fq.gz ..
    >>>

    output {
        File merged_fastq = "~{prefix}.fq.gz"
    }

    #########################
    Int memory = 8
    Int disk_size = 3 * ceil(size(fastqs, "GiB"))

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             memory,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
