version 1.0

import "../../structs/Structs.wdl"

workflow RunPBSV {

    meta {
        description: "Run PBSV to call SVs from a BAM file."
    }

    parameter_meta {
        bam:               "input BAM from which to call SVs"
        bai:               "index accompanying the BAM"
        is_ccs:            "if input BAM is CCS reads"
        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"
        prefix:            "prefix for output"
        zones:             "zones to run in"
        tandem_repeat_bed: "BED file containing TRF finder results (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"
    }

    input {
        File bam
        File bai
        Boolean is_ccs

        File ref_fasta
        File ref_fasta_fai
        String prefix

        Array[String] zones

        File? tandem_repeat_bed
    }

    call Discover {
        input:
            bam               = bam,
            bai               = bai,
            ref_fasta         = ref_fasta,
            ref_fasta_fai     = ref_fasta_fai,
            tandem_repeat_bed = tandem_repeat_bed,
            prefix            = prefix,
            zones             = zones
    }

    call Call {
        input:
            svsigs        = [ Discover.svsig ],
            ref_fasta     = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ccs           = is_ccs,
            prefix        = prefix,
            zones         = zones
    }

    output {
        File vcf = Call.vcf
    }
}

task Discover {
    input {
        File bam
        File bai
        File ref_fasta
        File ref_fasta_fai
        File? tandem_repeat_bed
        String? chr
        String prefix
        Array[String] zones
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

    Int MINIMAL_DISK = 500
    Boolean is_big_bam = size(bam, "GB") > 100
    Int inflation_factor = if (is_big_bam) then 5 else 2
    Int disk_size = inflation_factor * (ceil(size([bam, bai, ref_fasta, ref_fasta_fai], "GB")) + 1)
    Int runtime_disk_size = if disk_size < MINIMAL_DISK then MINIMAL_DISK else disk_size

    String fileoutput = if defined(chr) then "~{prefix}.~{chr}.svsig.gz" else "~{prefix}.svsig.gz"

    command <<<
        set -euxo pipefail

        pbsv discover \
            ~{if defined(tandem_repeat_bed) then "--tandem-repeats ~{tandem_repeat_bed}" else ""} \
            ~{bam} \
            ~{fileoutput}
    >>>

    output {
        File svsig = "~{fileoutput}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          if(defined(chr)) then 8 else 32,
        mem_gb:             if(defined(chr)) then 32 else 128,
        disk_gb:            runtime_disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sv:0.1.8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        zones: zones
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
        Boolean ccs
        String prefix
        Array[String] zones
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

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        pbsv call -j $num_core --log-level INFO ~{true='--ccs' false='' ccs} \
            ~{ref_fasta} \
            ~{sep=' ' svsigs} \
            ~{prefix}.pbsv.pre.vcf

        cat ~{prefix}.pbsv.pre.vcf | grep -v -e '##fileDate' > ~{prefix}.pbsv.vcf
    >>>

    output {
        File vcf = "~{prefix}.pbsv.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             96,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sv:0.1.8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        zones: zones
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

