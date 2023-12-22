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
        output_bucket:     "cloud path for output files"
        tandem_repeat_bed: "BED file containing TRF finder results (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"
        regions:           "genomic regions to process"
        n_tasks:           "number of tasks into which a whole genome is divided"
    }

    input {
        File bam
        File bai
        Boolean is_ccs

        File ref_fasta
        File ref_fasta_fai
        String prefix
        String? output_bucket

        File? tandem_repeat_bed
        Array[String]? regions
        Int n_tasks = 1
    }

    call PBSV {
        input:
            bam               = bam,
            bai               = bai,
            is_ccs            = is_ccs,
            ref_fasta         = ref_fasta,
            ref_fasta_fai     = ref_fasta_fai,
            prefix            = prefix,
            output_bucket     = output_bucket,
            tandem_repeat_bed = tandem_repeat_bed,
            regions           = regions,
            n_tasks           = n_tasks
    }

    output {
        File vcf = PBSV.vcf
        File tbi = PBSV.tbi
    }
}

task PBSV {
    input {
        File bam
        File bai
        Boolean is_ccs
        File ref_fasta
        File ref_fasta_fai
        String prefix
        String? output_bucket
        File? tandem_repeat_bed
        Array[String]? regions
        Int n_tasks
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:               {description: "input BAM from which to call SVs",
                            localization_optional: true}
        bai:               {description: "index accompanying the BAM",
                            localization_optional: true}
        is_ccs:            "if input BAM is CCS reads"
        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"
        prefix:            "prefix for output"
        output_bucket:     "cloud path for output files"
        tandem_repeat_bed: "BED file containing TRF finder (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"
        regions:           "genomic regions to process"
        n_tasks:           "number of tasks into which a whole genome is divided"
    }

    Int MINIMAL_DISK = 500
    Float bam_size = size(bam, "GB")/n_tasks
    Int inflation_factor = if (bam_size > 100) then 5 else 2
    Int disk_size = inflation_factor * (ceil(bam_size + size(ref_fasta, "GB")) + 1)
    Int runtime_disk_size = if disk_size < MINIMAL_DISK then MINIMAL_DISK else disk_size

    command <<<
        set -euxo pipefail

        if ~{defined(regions)}; then
            GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
            export GCS_OAUTH_TOKEN
            samtools view -h1X --write-index -o "local.bam##idx##local.bam.bai" "~{bam}" "~{bai}" "~{sep='" "' regions}"
        else
            gcloud storage cp "~{bam}" local.bam
        fi

        env # testing: having trouble finding pbsv

        pbsv discover \
            ~{"--tandem-repeats " + tandem_repeat_bed} local.bam local.svsig.gz

        pbsv call -j "$(nproc)" --log-level INFO ~{true='--ccs' false='' is_ccs} \
            "~{ref_fasta}" \
            local.svsig.gz \
            "~{prefix}.pbsv.vcf"

        vcfName="~{prefix}.pbsv.vcf.gz"
        grep -v '##fileDate' < "~{prefix}.pbsv.vcf" | bgzip > "$vcfName"
        tabix -p vcf "$vcfName"

        if ~{defined(output_bucket)}; then
            outDir=$(echo "~{output_bucket}" | sed 's+/?$+/+')
            gcloud storage cp "$vcfName" "${vcfName}.tbi" "$outDir"
        fi
    >>>

    output {
        File vcf = "~{prefix}.pbsv.vcf.gz"
        File tbi = "~{prefix}.pbsv.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          if(defined(regions)) then 8 else 32,
        mem_gb:             if(defined(regions)) then 64 else 128,
        disk_gb:            runtime_disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pbsv:0.1.1"
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
