version 1.0

import "../../structs/Structs.wdl"

task AnnotateAdapters {
    input {
        File bam
        Int read_end_length = 500
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 2
    Int disk_size = 4*ceil(size(bam, "GB"))

    String output_name = basename(bam, ".bam")

    command <<<
        set -euxo pipefail

        python3 /lrma/tool.py \
            --bam=~{bam} \
            --adapter=/lrma/adapter_sequence.fasta \
            --reverse-adapter=/lrma/reverse_adapter_sequence.fasta \
            --whitelist-10x=/lrma/3M-february-2018.txt \
            --name=~{output_name}_annotated \
            --read-end-length=~{read_end_length} \
            --record-umis \
            --ssw-path /lrma/ssw/ \
            --starcode-path /lrma/starcode-master/starcode

        /opt/conda/envs/10x_tool/bin/samtools fastq -T ZA,CR,ZU,CB ~{output_name}_annotated.bam | gzip > ~{output_name}_annotated.fastq.gz
    >>>

    output {
        File annotated_bam  = "~{output_name}_annotated.bam"
        File annotated_fq   = "~{output_name}_annotated.fastq.gz"
        File barcode_stats  = "~{output_name}_annotated_barcode_stats.tsv"
        File starcode_stats = "~{output_name}_annotated_starcode.tsv"
        File stats          = "~{output_name}_annotated_stats.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-10x:0.1.9"
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
