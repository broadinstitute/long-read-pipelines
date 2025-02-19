version 1.0

import "../../structs/Structs.wdl"

workflow RunPBSV {

    meta {
        description: "Run PBSV to call SVs from a BAM file."
    }

    parameter_meta {
        bam:               "input BAM from which to call SVs"
        bai:               "index accompanying the BAM"
        is_hifi:           "if input BAM is HiFi reads"
        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"
        prefix:            "prefix for output"
        zones:             "zones to run in"
        tandem_repeat_bed: "BED file containing TRF finder results (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"
    }

    input {
        File bam
        File bai
        Boolean is_hifi

        File ref_fasta
        File ref_fasta_fai
        String prefix

        String zones

        File? tandem_repeat_bed
        Boolean is_ont = false
    }

    call Discover {
        input:
            bam               = bam,
            bai               = bai,
            is_hifi           = is_hifi,
            ref_fasta         = ref_fasta,
            ref_fasta_fai     = ref_fasta_fai,
            tandem_repeat_bed = tandem_repeat_bed,
            prefix            = prefix,
            zones             = zones,
            is_ont            = is_ont
    }

    call Call {
        input:
            svsigs        = [ Discover.svsig ],
            ref_fasta     = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            is_hifi       = is_hifi,
            prefix        = prefix,
            zones         = zones,
            is_ont        = is_ont
    }

    output {
        File vcf = Call.vcf
        File tbi = Call.tbi
    }
}

task Discover {
    input {
        File bam
        File bai
        Boolean is_hifi
        File ref_fasta
        File ref_fasta_fai
        File? tandem_repeat_bed
        String? chr
        String prefix
        String zones
        Boolean is_ont = false
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:               "input BAM from which to call SVs"
        bai:               "index accompanying the BAM"
        is_hifi:           "if input BAM is HiFi reads"
        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"
        tandem_repeat_bed: "BED file containing TRF finder (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"
        chr:               "chr on which to call variants"
        prefix:            "prefix for output"
    }

    String fileoutput = if defined(chr) then "~{prefix}.~{chr}.svsig.gz" else "~{prefix}.svsig.gz"

    command <<<
        set -euxo pipefail

        # pbsv, for ONT inputs, could fail with a really strange error that looks like the following
        # >|> 20230828 01:52:02.907 -|- FATAL -|- Run -|- 0x7f704219e4c0|| -|- pbsv discover ERROR: map::at
        # upon investigation, it's most likely caused by fields in the RG lines in the header that are not ID, SM, PU
        # so we fix it here
        # hastag facts-in-life
        if ~{is_ont}; then
            samtools view --no-PG -H ~{bam} > orig.header.txt
            grep -v "^@SQ" orig.header.txt

            # keep the bare minimum info in @RG
            grep "^@RG" orig.header.txt > RG.lines
            cat RG.lines
            while IFS= read -r line; do
                id_field=$(echo "${line}" | tr '\t' '\n' | grep "^ID:")
                pu_field=$(echo "${line}" | tr '\t' '\n' | grep "^PU:")
                sm_field=$(echo "${line}" | tr '\t' '\n' | grep "^SM:")
                echo -e "@RG\t${id_field}\t${pu_field}\t${sm_field}" >> fixed.rg.lines
            done < RG.lines
            cat fixed.rg.lines

            # patch back the header
            # this is to follow the order we observe typically in BAM headers: @HD first, @SQ lines, @RG lines, then the rest (usually @PG lines)
            grep -v "^@RG" orig.header.txt > non.RG.lines
            cat <(head -n1 non.RG.lines) \
                <(grep "^@SQ" non.RG.lines) \
                fixed.rg.lines \
                <(tail +2 non.RG.lines | grep -v "^@HD" | grep -v "^@SQ") \
            > fixed.header.txt
            cat fixed.header.txt

            date
            samtools reheader --no-PG fixed.header.txt ~{bam} | samtools view -@1 --no-PG -o tmp.bam
            date

            samtools view -H tmp.bam | grep "^@RG"
            mv tmp.bam ~{bam}
        fi
        pbsv discover \
            ~{true='--hifi' false='' is_hifi} \
            ~{if defined(tandem_repeat_bed) then "--tandem-repeats ~{tandem_repeat_bed}" else ""} \
            ~{bam} \
            ~{fileoutput}
    >>>

    output {
        File svsig = "~{fileoutput}"
    }

    #########################

    Int MINIMAL_DISK = 20
    Boolean is_big_bam = size(bam, "GB") > 100
    Int inflation_factor = if (is_big_bam) then 3 else 2
    Int disk_size = inflation_factor * (ceil(size([bam, bai, ref_fasta, ref_fasta_fai], "GB")) + 1)
    Int runtime_disk_size = if disk_size < MINIMAL_DISK then MINIMAL_DISK else disk_size

    String disk_type = if is_ont then "SSD" else "HDD"

    Int num_cores = if(defined(chr)) then 4 else 32
    Int memory = 2 * num_cores

    RuntimeAttr default_attr = object {
        cpu_cores:          num_cores,
        mem_gb:             memory,
        disk_gb:            runtime_disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-smrttools:12.0.0.176214"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        zones: zones
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
        Boolean is_hifi
        String prefix
        String zones
        Boolean DEBUG = false
        Boolean is_ont = false
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        svsigs:            "per-chromosome *.svsig.gz files"
        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"
        is_hifi:           "if input BAM is HiFi reads"
        prefix:            "prefix for output"
    }


    command <<<
        set -euxo pipefail

        pbsv call \
            -j 0 \
            --log-level ~{true='INFO' false='WARN' DEBUG} \
            --log-file pbsv.call.log \
            ~{true='--hifi' false='' is_hifi} \
            ~{ref_fasta} \
            ~{sep=' ' svsigs} \
            ~{prefix}.pbsv.pre.vcf

        # some trivial postprocessing
        cat ~{prefix}.pbsv.pre.vcf | grep -v -e '##fileDate' > ~{prefix}.pbsv.vcf
        bgzip -c ~{prefix}.pbsv.vcf > ~{prefix}.pbsv.vcf.gz
        tabix -p vcf ~{prefix}.pbsv.vcf.gz
    >>>

    output {
        File call_log = "pbsv.call.log"  # make sure this is always the top, so that in case something goes wrong, we still get the log de-localized
        File vcf = "~{prefix}.pbsv.vcf.gz"
        File tbi = "~{prefix}.pbsv.vcf.gz.tbi"
    }

    #########################
    Int disk_size = 2*ceil(size(svsigs, "GiB") + size([ref_fasta, ref_fasta_fai], "GiB"))

    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             128,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-smrttools:12.0.0.176214"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        zones: zones
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
