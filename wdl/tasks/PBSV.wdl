version 1.0

##########################################################################################
# This pipeline calls SVs on an input LR BAM using various known SV algorithms
# that are specifically designed to work with long read data.
# Each individual task/algo. is directly callable, if so desired.
##########################################################################################

import "Structs.wdl"

task Discover {
    input {
        String bam
        String bai

        File ref_fasta
        File ref_fasta_fai
        File? tandem_repeat_bed

        String chr
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:               "input BAM from which to call SVs"
        bai:               "index accompanying the BAM"

        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"
        tandem_repeat_bed: "BED file containing TRF finder (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"

        chr:               "chr on which to call variants"
        prefix:            "prefix for output"
    }

    Int disk_size = 2*ceil(size([bam, bai, ref_fasta, ref_fasta_fai, tandem_repeat_bed], "GB"))

    # purely experiential
    Int memory = if (ceil(size(bam, "GiB")) > 20) then 96 else 64
    Int cpus = ceil( memory / 6 ) # a range of, approximately, [1,6] ratio between mem/cpu allowed from cloud service provider

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -hb ~{bam} ~{chr} > ~{chr}.bam
        samtools index ~{chr}.bam

        pbsv discover \
            ~{if defined(tandem_repeat_bed) then "--tandem-repeats ~{tandem_repeat_bed}" else ""} \
            ~{chr}.bam \
            ~{prefix}.~{chr}.svsig.gz
    >>>

    output {
        File svsig = "~{prefix}.~{chr}.svsig.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sv:0.1.8"
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

task Call {
    input {
        Array[File] svsigs

        File ref_fasta
        File ref_fasta_fai

        Boolean ccs = false

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        svsigs:            "per-chromosome *.svsig.gz files"

        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"

        ccs:               "use optimizations for CCS data"
        prefix:            "prefix for output"
    }

    Int disk_size = 2*ceil(size(svsigs, "GiB") + size([ref_fasta, ref_fasta_fai], "GiB"))

    command <<<
        set -euxo pipefail

        pbsv call ~{true='--ccs' false='' ccs} ~{ref_fasta} ~{sep=' ' svsigs} ~{prefix}.pbsv.pre.vcf

        cat ~{prefix}.pbsv.pre.vcf | grep -v -e '##fileDate' > ~{prefix}.pbsv.vcf
    >>>

    output {
        File vcf = "~{prefix}.pbsv.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sv:0.1.8"
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
