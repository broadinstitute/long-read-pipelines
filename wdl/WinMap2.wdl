version 1.0

import "tasks/Structs.wdl"

workflow WinMap2 {
    input {
        Array[File] queries
        File target_fasta

        String map_preset

        String sam_out_prefix
    }

    parameter_meta {
        queries:        "query sequences to be mapped and aligned"
        target_fasta:   "target fasta to map against"
        map_preset:     "controled vocabulary (ONT, PB, ASM) for mapping long reads or assemblies"
        sam_out_prefix: "prefix for output BAM"
    }

    call MerylPrepRef {
        input:
        target_fasta = target_fasta,
        kmer_size = if (map_preset=="ASM") then 19 else 15  # straight from README of winnowmap
    }

    if (length(queries) == 1) {
        call Winnowmap_map {
            input:
            query = select_first(queries),
            target_fasta = target_fasta,
            target_meryl_repetitive_info = MerylPrepRef.repetitive_txt,
            map_preset = map_preset,
            sort_and_index = true,
            sam_out_prefix = sam_out_prefix
        }
    }

    if (length(queries) > 1) {
        scatter (q in queries) {
            call Winnowmap_map as map {
                input:
                query = q,
                target_fasta = target_fasta,
                target_meryl_repetitive_info = MerylPrepRef.repetitive_txt,
                map_preset = map_preset,
                sort_and_index = false,
                sam_out_prefix = sam_out_prefix
            }
        }

        call MergeSortAndIndex as merged {
            input: unsorted_bams = map.aligned_bam, sam_out_prefix = sam_out_prefix
        }
    }

    output {
        File winnowmap_bam = select_first([Winnowmap_map.aligned_bam, merged.aligned_bam])
        File winnowmap_bai = select_first([Winnowmap_map.aligned_bai, merged.aligned_bai])
    }
}

task MerylPrepRef {
    # for preping a reference when using winnowmap
    # todo: we can think about optionally running this for regularly used references

    input {
        File target_fasta
        Int kmer_size
        Float distinct = 0.998

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        target_fasta: "reference fasta"
        kmer_size:    "argument for k for meryl count"
        distinct:     "argument for 'distinct' for meryl print greater-than"
    }

    command <<<
        set -euxo pipefail

        meryl count k="~{kmer_size}" output merylDB "~{target_fasta}"
        meryl print greater-than distinct="~{distinct}" merylDB > "repetitive_k~{kmer_size}.txt"
    >>>

    output {
        File repetitive_txt = "repetitive_k~{kmer_size}.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             10,
        disk_gb:            1 + 2*ceil(size(target_fasta, "GB")),
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-winnowmap:2.03"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    runtime_attr.cpu_cores
        memory:                 runtime_attr.mem_gb + " GiB"
        disks: "local-disk " +  runtime_attr.disk_gb + " HDD"
        bootDiskSizeGb:         runtime_attr.boot_disk_gb
        preemptible:            runtime_attr.preemptible_tries
        maxRetries:             runtime_attr.max_retries
        docker:                 "~{runtime_attr.docker}"
    }
}

task Winnowmap_map {
    input {
        File query

        File target_fasta
        File target_meryl_repetitive_info

        String map_preset
        Boolean sort_and_index

        String sam_out_prefix
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        query:                        "query sequences to be mapped and aligned"
        target_fasta:                 "target fasta to map against"
        target_meryl_repetitive_info: "meryl repetitive infomation file"
        sort_and_index:               "whether to sort & index the SAM file"
        sam_out_prefix:               "prefix for output BAM"
    }

    Int prefered_threads = 4

    command <<<
        set -euxo pipefail

        if [[ "~{map_preset}" == 'ASM' ]]; then 
            mapping_params='asm20'
        elif [[ "~{map_preset}" == 'ONT' ]]; then
            mapping_params='map-ont'
        elif [[ "~{map_preset}" == 'PB' ]]; then
            mapping_params='map-pb'
        else
            echo "mapping preset need to be one of (ONT, PB, ASM). Received ~{map_preset}"
            exit 1
        fi

        if [[ "~{sort_and_index}" ]]; then
            winnowmap \
                -W "~{target_meryl_repetitive_info}" \
                -t ~{prefered_threads} \
                -ax "${mapping_params}" \
                "~{target_fasta}" \
                "~{query}" | \
                samtools sort -o ~{sam_out_prefix}.bam -@ ~{prefered_threads} -
            samtools index -@ ~{prefered_threads} ~{sam_out_prefix}.bam
        else
            winnowmap \
                -W "~{target_meryl_repetitive_info}" \
                -t ~{prefered_threads} \
                -ax "${mapping_params}" \
                "~{target_fasta}" \
                "~{query}" | \
                samtools view -o ~{sam_out_prefix}.bam -@ ~{prefered_threads} -
        fi
    >>>

    output {
        File  aligned_bam = "~{sam_out_prefix}.bam"
        File? aligned_bai = "~{sam_out_prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          prefered_threads + 2,
        mem_gb:             30,
        disk_gb:            1 + 3*ceil(size(query, "GB") + size(target_fasta, "GB")),
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-winnowmap:2.03"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    runtime_attr.cpu_cores
        memory:                 runtime_attr.mem_gb + " GiB"
        disks: "local-disk " +  runtime_attr.disk_gb + " HDD"
        bootDiskSizeGb:         runtime_attr.boot_disk_gb
        preemptible:            runtime_attr.preemptible_tries
        maxRetries:             runtime_attr.max_retries
        docker:                 "~{runtime_attr.docker}"
    }
}

task MergeSortAndIndex {
    input {
        Array[File] unsorted_bams
        String sam_out_prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        samtools merge -o unsorted.bam ~{sep=' ' unsorted_bams}
        rm -rf ~{sep=' ' unsorted_bams}
        samtools sort -o "~{sam_out_prefix}.bam" unsorted.bam
        samtools index "~{sam_out_prefix}.bam"
    >>>

    output {
        File aligned_bam = "~{sam_out_prefix}.bam"
        File aligned_bai = "~{sam_out_prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             30,
        disk_gb:            1 + 2*size(unsorted_bams, "GB"),
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-winnowmap:2.03"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    runtime_attr.cpu_cores
        memory:                 runtime_attr.mem_gb + " GiB"
        disks: "local-disk " +  runtime_attr.disk_gb + " HDD"
        bootDiskSizeGb:         runtime_attr.boot_disk_gb
        preemptible:            runtime_attr.preemptible_tries
        maxRetries:             runtime_attr.max_retries
        docker:                 "~{runtime_attr.docker}"
    }
}
