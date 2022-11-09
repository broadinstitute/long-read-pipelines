version 1.0

import "Structs.wdl"

task MergePairedEndReads {
    input {
        File illumina_fq1
        File illumina_fq2

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            25,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "quay.io/biocontainers/pear:0.9.6--h67092d7_8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int num_threads = round(select_first([runtime_attr.cpu_cores, default_attr.cpu_cores]))

    command <<<
        pear -j ~{num_threads} -f ~{illumina_fq1} -r ~{illumina_fq2} -o pear
    >>>

    output {
        File merged = "pear.assembled.fastq"
        File split1 = "pear.unassembled.forward.fastq"
        File split2 = "pear.unassembled.reverse.fastq"
        File discarded = "pear.discarded.fastq"
    }

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

task McCortexBuild {
    input {
        String sample_id
        Int k

        File illumina_fq1
        File? illumina_fq2

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            25,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "quay.io/biocontainers/mccortex:1.0--hd03093a_5"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int max_mem = round(select_first([runtime_attr.mem_gb, default_attr.mem_gb])) - 2  # 2 GB buffer
    Int num_threads = round(select_first([runtime_attr.cpu_cores, default_attr.cpu_cores]))

    String flag_str = if defined(illumina_fq2) && illumina_fq2 != "" then "--seq2 " else "--seq "
    String flag_sep = if defined(illumina_fq2) && illumina_fq2 != "" then ":" else ""

    command <<<
        mccortex ~{k} build -t ~{num_threads} -m ~{max_mem}G -k ~{k} --sample ~{sample_id} \
            ~{flag_str}~{illumina_fq1}~{flag_sep}~{illumina_fq2} \
            ~{sample_id}.ctx
        mccortex ~{k} clean -t ~{num_threads} -m ~{max_mem}G --out ~{sample_id}.cleaned.ctx ~{sample_id}.ctx

        mccortex ~{k} inferedges -t ~{num_threads} -m ~{max_mem}G ~{sample_id}.ctx
        mccortex ~{k} inferedges -t ~{num_threads} -m ~{max_mem}G ~{sample_id}.cleaned.ctx
    >>>

    output {
        File graph = "~{sample_id}.ctx"
        File graph_cleaned = "~{sample_id}.cleaned.ctx"
    }

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

task McCortexLinksForRef {
    input {
        Int k

        File mccortex_graph

        String ref_id
        File ref_fasta

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            25,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "quay.io/biocontainers/mccortex:1.0--hd03093a_5"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int max_mem = round(select_first([runtime_attr.mem_gb, default_attr.mem_gb])) - 2  # 2 GB buffer
    Int num_threads = round(select_first([runtime_attr.cpu_cores, default_attr.cpu_cores]))

    command <<<
        mccortex ~{k} thread -t ~{num_threads} -m ~{max_mem} --seq ~{ref_fasta} --out ~{ref_id}.ctp.gz ~{mccortex_graph}
    >>>

    output {
        File mccortex_links = "~{ref_id}.ctp.gz"
    }

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

task McCortexLinksForReads {
    input {
        String sample_id
        Int k

        File merged_fq
        File illumina_fq1
        File illumina_fq2

        File mccortex_graph

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            25,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "quay.io/biocontainers/mccortex:1.0--hd03093a_5"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int max_mem = round(select_first([runtime_attr.mem_gb, default_attr.mem_gb])) - 2  # 2 GB buffer
    Int num_threads = round(select_first([runtime_attr.cpu_cores, default_attr.cpu_cores]))

    command <<<
        # Links from merged paired-end reads and then from non-merged reads
        mccortex ~{k} thread -t ~{num_threads} -m ~{max_mem} --seq ~{merged_fq} --out ~{sample_id}.merged.raw.ctp.gz ~{mccortex_graph}
        mccortex ~{k} thread -t ~{num_threads} -m ~{max_mem} --seq2 ~{illumina_fq1}:~{illumina_fq2} --out ~{sample_id}.other.raw.ctp.gz ~{mccortex_graph}

        # Combine them in a single link file
        mccortex ~{k} pjoin ~{sample_id}.merged.raw.ctp.gz --out ~{sample_id}.raw.ctp.gz ~{sample_id}.other.raw.ctp.gz

        # Prune low coverage links
        mccortex ~{k} links -T link.stats.txt -L 1000 ~{sample_id}.raw.ctp.gz
        LINK_THRESH=`grep 'suggested_cutoff=' link.stats.txt | grep -oE '[0-9,]+$'`
        mccortex ~{k} --clean $LINK_THRESH --out ~{sample_id}.ctp.gz ~{sample_id}.raw.ctp.gz
    >>>

    output {
        File mccortex_links = "~{sample_id}.ctp.gz"
    }

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

task McCortexAssemble {
    input {
        String sample_id
        Int k

        File mccortex_graph
        Array[File] ref_links
        File sample_links

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            25,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "quay.io/biocontainers/mccortex:1.0--hd03093a_5"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int max_mem = round(select_first([runtime_attr.mem_gb, default_attr.mem_gb])) - 2  # 2 GB buffer
    Int num_threads = round(select_first([runtime_attr.cpu_cores, default_attr.cpu_cores]))

    command <<<
        # Join ref links and sample links
        mccortex ~{k} pjoin $(< ~{write_lines(ref_links)}) -o all_links.ctp.gz ~{sample_links}
        mccortex ~{k} popbubbles -t ~{num_threads} -m ~{max_mem} --out popped_bubbles.ctx ~{mccortex_graph}
        mccortex ~{k} contigs -t ~{num_threads} -m ~{max_mem} -p all_links.ctp.gz -o ~{sample_id}.contigs.fasta popped_bubbles.ctx
    >>>

    output {
        File contigs = "~{sample_id}.contigs.fasta"
    }

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
