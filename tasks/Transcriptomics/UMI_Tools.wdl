version 1.0

import "../Structs.wdl"

task Run_Group {

    meta {
        description : "Run umi-tools group on a bam file."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File aligned_transcriptome_reads
        File aligned_transcriptome_reads_index

        String gene_tag = "XG"
        String cell_barcode_tag = "CB"
        String umi_tag = "ZU"

        Boolean do_per_cell = true

        String prefix = "umi_tools_group"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 20*ceil(size(aligned_transcriptome_reads, "GB") + size(aligned_transcriptome_reads_index, "GB"))

    String per_cell_args = if do_per_cell then " --per-cell --cell-tag " + cell_barcode_tag + " " else ""

    String memory_log_file = "memory_use.txt"

    command <<<

        # Set up memory logging daemon:
        MEM_LOG_INTERVAL_s=5
        DO_MEMORY_LOG=true
        while $DO_MEMORY_LOG ; do
            date
            date +%s
            cat /proc/meminfo
            sleep $MEM_LOG_INTERVAL_s
        done >> ~{memory_log_file} &
        mem_pid=$!

        set -euxo pipefail

        # Run umi-tools group:
        umi_tools group \
          --buffer-whole-contig \
          --no-sort-output \
          --per-gene \
          ~{per_cell_args} \
          --gene-tag ~{gene_tag} \
          --extract-umi-method tag \
          --umi-tag ~{umi_tag} \
          -I ~{aligned_transcriptome_reads} \
          --group-out=~{prefix}.tsv \
          --output-bam \
          --log=~{prefix}.log > ~{prefix}.bam


        # Stop the memory daemon softly.  Then stop it hard if it's not cooperating:
        set +e
        DO_MEMORY_LOG=false
        sleep $(($MEM_LOG_INTERVAL_s  * 2))
        kill -0 $mem_pid &> /dev/null
        if [ $? -ne 0 ] ; then
            kill -9 $mem_pid
        fi
    >>>

    output {
        File output_bam = "~{prefix}.bam"
        File output_tsv = "~{prefix}.tsv"
        File memory_log = "~{memory_log_file}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             64,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
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