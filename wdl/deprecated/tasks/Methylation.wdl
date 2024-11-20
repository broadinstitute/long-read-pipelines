version 1.0

##########################################################################################
# Workflow that runs F5C rewrite of Nanopolish to call methylation.
##########################################################################################

import "../../structs/Structs.wdl"

workflow Methylation {
    input {
        Array[File] fast5s
        Array[File] fastqs
        File? sequencing_summary

        File bam
        File bai

        File ref_fasta
        File ref_fai

        String prefix = "out"
    }

    parameter_meta {
        fast5s: "input raw Nanopore signal data"
        fastqs: "basecalled fastqs"
        sequencing_summary: "summary file from Nanopore basecaller"

        bam: "aligned reads"
        bai: "aligned reads index"

        ref_fasta: "reference to which the reads were aligned"
        ref_fai:   "index accompanying the reference"

        prefix: "[default-valued] prefix for output files"
    }

    call CombineFastqs { input: fastqs = fastqs }

    call CallMethylation {
        input:
            fast5s = fast5s,
            sequencing_summary = sequencing_summary,

            fastq = CombineFastqs.fastq,

            bam = bam,
            bai = bai,

            ref_fasta = ref_fasta,
            ref_fai = ref_fai,

            prefix = prefix
    }

    output {
        File index = CallMethylation.index
        File index_fai = CallMethylation.index_fai
        File index_gzi = CallMethylation.index_gzi
        File index_readdb = CallMethylation.index_readdb

        File results_tsv = CallMethylation.results_tsv
        File freq_tsv = CallMethylation.freq_tsv
    }
}

task CombineFastqs {
    input {
        Array[File] fastqs
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        fastqs: "basecalled fastq files"
        prefix: "[default-valued] prefix for combined fastq"
    }

    Int disk_size = 1 + 3*ceil(size(fastqs, "GB"))

    command <<<
        set -euxo pipefail

        find . -type f -name '*.fastq' -exec cat {} + > ~{prefix}.fastq
    >>>

    output {
        File fastq = "~{prefix}.fastq"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
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

task CallMethylation {
    input {
        Array[File] fast5s
        File? sequencing_summary

        File fastq

        File bam
        File bai

        File ref_fasta
        File ref_fai

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        fast5s: "input raw Nanopore signal data"
        sequencing_summary: "summary file from Nanopore basecaller"

        fastq: "basecalled fastqs concatenated into a single fastq file"

        bam: "aligned reads"
        bai: "aligned reads index"

        ref_fasta: "reference to which the reads were aligned"
        ref_fai:   "index accompanying the reference"

        prefix: "[default-valued] prefix for output files"
    }

    Int disk_size = 2*ceil(size(fast5s, "GB") + size(select_all([sequencing_summary, fastq, bam, bai, ref_fasta, ref_fai]), "GB"))

    command <<<
        set -euxo pipefail

        lspci | grep -i "vga\|3d\|display"
        lspci -nnk | grep -iA2 "vga\|3d\|display"
        nvcc --version

        FAST5_DIR=$(dirname ~{fast5s[0]})

        f5c index ~{if defined(sequencing_summary) then "-s ~{sequencing_summary}" else "" } -d $FAST5_DIR ~{fastq}
        f5c call-methylation -b ~{bam} -g ~{ref_fasta} -r ~{fastq} --iop 32 -K 2048 -B 9.3M > ~{prefix}.results.tsv
        f5c meth-freq -i ~{prefix}.results.tsv > ~{prefix}.freq.tsv
    >>>

    output {
        File index = "~{fastq}.index"
        File index_fai = "~{fastq}.index.fai"
        File index_gzi = "~{fastq}.index.gzi"
        File index_readdb = "~{fastq}.index.readdb"

        File results_tsv = "~{prefix}.results.tsv"
        File freq_tsv = "~{prefix}.freq.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-methylation:0.1.2"
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
        gpuType:                "nvidia-tesla-p100"
        gpuCount:               1
        nvidiaDriverVersion:    "418.152.00"
        zones:                  ["us-east1-b", "us-east1-c"]
        #cpuPlatform:            "Intel Skylake"
    }
}

task FreqMerge {
    input {
        Array[File] freq_tsvs
        String prefix = "combined"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        freq_tsvs: "methylation frequency tsv files"
        prefix: "[default-valued] prefix for combined methylation frequency tsv"
    }

    Int disk_size = 1 + 2*ceil(size(freq_tsvs, "GB"))

    command <<<
        set -euxo pipefail

        f5c freq-merge -o ~{prefix}.tsv -n ~{length(freq_tsvs)} -f ~{sep=' ' freq_tsvs}
    >>>

    output {
        File freq_tsv = "~{prefix}.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-methylation:0.1.0"
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
