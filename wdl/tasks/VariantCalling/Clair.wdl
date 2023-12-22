version 1.0

#######################################################
# This pipeline calls small variants using DeepVariant.
#######################################################

import "../../structs/Structs.wdl"

task Clair {

    meta {
        description: "Call variants using Clair3."
    }

    parameter_meta {
        bam:             {description: "input BAM from which to call variants",
                          localization_optional: true}
        bai:             {description: "index accompanying the BAM",
                          localization_optional: true}

        ref_fasta:       "reference to which the BAM was aligned"
        ref_fasta_fai:   "index accompanying the reference"

        sites_vcf:       "sites VCF"
        sites_vcf_tbi:   "sites VCF index"

        regions:         "genomic regions on which to call variants"
        preset:          "calling preset (CCS, ONT)"

        runtime_attr_override: "override the default runtime attributes"
    }

    input {
        File bam
        File bai

        String prefix
        String? output_bucket

        File ref_fasta
        File ref_fasta_fai

        File? sites_vcf
        File? sites_vcf_tbi

        Array[String] regions
        String preset

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(select_all([bam, bai, ref_fasta, ref_fasta_fai, sites_vcf]), "GB"))
    String platform = if preset == "CCS" then "hifi" else "ont"

    command <<<
        set -euxo pipefail

        GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        export GCS_OAUTH_TOKEN

        samtools view -h1X --write-index -o "local.bam##idx##local.bam.bai" "~{bam}" "~{bai}" "~{sep='" "' regions}"

        SAMPLE=$(samtools view -H local.bam | sed '/^@RG/!d;s/.*	SM:\([^	]*\).*/\1/' | sed '2,$d')

        /opt/bin/run_clair3.sh \
            ~{"--vcf_fn=" + sites_vcf} \
            --bam_fn=local.bam \
            --ref_fn="~{ref_fasta}" \
            --threads="$(nproc)" \
            --platform="~{platform}" \
            --model_path="/opt/models/~{platform}" \
            --sample_name="$SAMPLE" \
            --gvcf \
            --include_all_ctgs \
            --output="./"

        # for chrM, Clair3 creates a header only vcf, copy it to gVCF as-is
        if [[ ! -f merge_output.gvcf.gz ]]; then cp "merge_output.vcf.gz" "merge_output.gvcf.gz"; fi

        baseName="~{prefix}.clair"
        vcfName="${baseName}.vcf.gz"
        tbiName="${baseName}.vcf.gz.tbi"
        gvcfName="${baseName}.gvcf.gz"
        gtbiName="${baseName}.gvcf.gz.tbi"
        mv merge_output.vcf.gz "$vcfName"
        mv merge_output.vcf.gz.tbi "$tbiName"
        mv merge_output.gvcf.gz "$gvcfName"
        mv merge_output.gvcf.gz.tbi "$gtbiName"
        if ~{defined(output_bucket)}; then
            outDir=$(echo "~{output_bucket}" | sed 's+/?$+/+')
            gcloud storage cp "$vcfName" "$tbiName" "$gvcfName" "$gtbiName" "$outDir"
        fi
    >>>

    output {
        File? pileup_vcf = "pileup.vcf.gz"
        File? pileup_vcf_tbi = "pileup.vcf.gz.tbi"
        File? full_alignment_vcf = "full_alignment.vcf.gz"
        File? full_alignment_tbi = "full_alignment.vcf.gz.tbi"

        # save both VCF and gVCF
        File vcf = "~{prefix}.clair.vcf.gz"
        File? vcf_tbi = "~{prefix}.clair.vcf.gz.tbi"
        File gvcf = "~{prefix}.clair.gvcf.gz"
        File? gvcf_tbi = "~{prefix}.clair.gvcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          36,
        mem_gb:             72,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-clair:2.1.2"
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
