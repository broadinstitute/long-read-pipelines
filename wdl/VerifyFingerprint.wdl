version 1.0

import "tasks/Structs.wdl"
import "tasks/Finalize.wdl" as FF

workflow VerifyFingerprint {

    meta {
        description: "A workflow to check correctness of metadata on a flowcell, by genotyping it's BAM generated with its metadata, against a 'truth' genotyped VCF."
    }

    input {
        File aligned_bam
        File aligned_bai
        String expt_type

        Int? artificial_baseQ_for_CLR = 10

        File fingerprint_vcf

        File ref_map_file

        String gcs_out_root_dir
    }

    parameter_meta {
        aligned_bam:        "GCS path to aligned BAM file, supposed to be of the same sample as from the fingerprinting VCF"
        expt_type:          "There will be special treatment for 'CLR' data (minimum base quality for bases used when computing a fingerprint)"
        artificial_baseQ_for_CLR: "An artificial value for CLR reads used for fingerprint verification (CLR reads come with all 0 base qual)"

        fingerprint_vcf:    "Single sample Fingerprint VCF file from local database"


        ref_map_file:       "table indicating reference sequence and auxillary file locations"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/VerifyFingerprint"

    call GetVCFSampleName {
        input:
            fingerprint_vcf = fingerprint_vcf
    }

    call FilterGenotypesVCF {
        input:
            fingerprint_vcf = fingerprint_vcf
    }

    if (expt_type!='CLR') {
        call CheckFingerprint {
            input:
                aligned_bam     = aligned_bam,
                aligned_bai     = aligned_bai,
                fingerprint_vcf = FilterGenotypesVCF.read_to_use_vcf,
                vcf_sample_name = GetVCFSampleName.sample_name,
                haplotype_map   = ref_map['haplotype_map']
        }
    }
    if (expt_type=='CLR') {
        call ExtractRelevantCLRReads {
            input:
                aligned_bam     = aligned_bam,
                aligned_bai     = aligned_bai,
                fingerprint_vcf = FilterGenotypesVCF.read_to_use_vcf,
        }
        call ResetCLRBaseQual {
            input:
                bam = ExtractRelevantCLRReads.relevant_reads,
                bai = ExtractRelevantCLRReads.relevant_reads_bai,
                arbitrary_bq = select_first([artificial_baseQ_for_CLR])
        }
        call CheckCLRFingerprint {
            input:
                aligned_bam     = ResetCLRBaseQual.barbequed_bam,
                aligned_bai     = ResetCLRBaseQual.barbequed_bai,
                min_base_q      = artificial_baseQ_for_CLR,
                fingerprint_vcf = FilterGenotypesVCF.read_to_use_vcf,
                vcf_sample_name = GetVCFSampleName.sample_name,
                haplotype_map   = ref_map['haplotype_map']
        }
    }

    File summary_metrics = select_first([CheckFingerprint.summary_metrics, CheckCLRFingerprint.summary_metrics])
    File detail_metrics  = select_first([CheckFingerprint.detail_metrics,  CheckCLRFingerprint.detail_metrics])

    call FF.FinalizeToFile as FinalizeFingerprintSummaryMetrics { input: outdir = outdir, file = summary_metrics }
    call FF.FinalizeToFile as FinalizeFingerprintDetailMetrics  { input: outdir = outdir, file = detail_metrics }

    Map[String, String] metrics_map = select_first([CheckFingerprint.metrics_map, CheckCLRFingerprint.metrics_map])

    output {
        Float lod_expected_sample = metrics_map['LOD_EXPECTED_SAMPLE']

        File fingerprint_metrics = FinalizeFingerprintSummaryMetrics.gcs_path
        File fingerprint_details = FinalizeFingerprintDetailMetrics.gcs_path
    }
}

task GetVCFSampleName {
    input {
        File fingerprint_vcf
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eux

        GREPCMD="grep"
        if [[ ~{fingerprint_vcf} =~ \.gz$ ]]; then
            GREPCMD="zgrep"
        fi
        "${GREPCMD}" \
            "^#CHROM" \
            ~{fingerprint_vcf} \
            | awk '{print $10}' \
            > sample_name.txt
    >>>

    output {
        String sample_name = read_string("sample_name.txt")
    }

    ###################
    runtime {
        cpu: 2
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"ubuntu:20.04"
    }
}

task FilterGenotypesVCF {
    input {
        File fingerprint_vcf
        Array[String] filters = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']
    }

    parameter_meta {
        filters: "An array of chromosome names to filter out when verifying fingerprints"
    }

    command <<<
        set -eux

        GREPCMD="grep"
        if [[ ~{fingerprint_vcf} =~ \.gz$ ]]; then
            GREPCMD="zgrep"
        fi
        "${GREPCMD}" \
            -v \
            -e ' placeholder ' \
            ~{true='-e' false='' length(filters) > 0} \
            ~{sep=" -e " filters} \
            ~{fingerprint_vcf}  \
            > fingerprint.fixed.vcf
    >>>

    output {
        File read_to_use_vcf = "fingerprint.fixed.vcf"
    }

    ###################
    runtime {
        cpu: 2
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"ubuntu:20.04"
    }
}

task CheckFingerprint {

    meta {
        description: "Uses Picard tool CheckFingerprint to verify if the samples in provided VCF and BAM arise from the same biological sample"
    }
    input {
        File aligned_bam
        File aligned_bai

        File fingerprint_vcf
        String vcf_sample_name

        File haplotype_map

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        aligned_bam:{
            description:  "GCS path to aligned BAM file, supposed to be of the same sample as from the fingerprinting VCF",
            localization_optional: true
        }

        fingerprint_vcf:    "Fingerprint VCF file from local database; note that sample name must be the same as in BAM"

        haplotype_map:      "table indicating reference sequence and auxillary file locations"
    }

    Int disk_size = ceil(size([fingerprint_vcf, haplotype_map], "GB"))
    String prefix = basename(aligned_bam, ".bam")

    command <<<
        set -eux

        gatk CheckFingerprint \
            --INPUT ~{aligned_bam} \
            --GENOTYPES ~{fingerprint_vcf} \
            --EXPECTED_SAMPLE_ALIAS ~{vcf_sample_name} \
            --HAPLOTYPE_MAP ~{haplotype_map} \
            --OUTPUT ~{prefix}

        grep -v '^#' ~{prefix}.fingerprinting_summary_metrics | \
            grep -A1 READ_GROUP | \
            awk '
                {
                    for (i=1; i<=NF; i++)  {
                        a[NR,i] = $i
                    }
                }
                NF>p { p = NF }
                END {
                    for(j=1; j<=p; j++) {
                        str=a[1,j]
                        for(i=2; i<=NR; i++){
                            str=str" "a[i,j];
                        }
                        print str
                    }
                }' | \
            sed 's/ /\t/' \
            > metrics_map.txt

        mv ~{prefix}.fingerprinting_summary_metrics \
            ~{prefix}.fingerprinting_summary_metrics.txt
        mv ~{prefix}.fingerprinting_detail_metrics \
            ~{prefix}.fingerprinting_detail_metrics.txt
    >>>

    output {
        File summary_metrics = "~{prefix}.fingerprinting_summary_metrics.txt"
        File detail_metrics = "~{prefix}.fingerprinting_detail_metrics.txt"
        Map[String, String] metrics_map = read_map("metrics_map.txt")
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:             2,
        mem_gb:                4,
        disk_gb:               disk_size,
        boot_disk_gb:          10,
        preemptible_tries:     3,
        max_retries:           2,
        docker:                "us.gcr.io/broad-gatk/gatk:4.2.0.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                   select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:        select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:           select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:            select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                select_first([runtime_attr.docker, default_attr.docker])
    }
}

task ExtractRelevantCLRReads {
    meta {
        description: "Based on genotyping (SNP) sites, extract reads that overlap those places"
    }
    input {
        File aligned_bam
        File aligned_bai

        File fingerprint_vcf
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        aligned_bam:{
            localization_optional: true
        }
    }

    command <<<

        set -eux

        grep -v "^#" ~{fingerprint_vcf} | \
            awk 'BEGIN {OFS="\t"} {print $1, $2-1, $2, $3}' \
            > genotyping.sites.bed

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        
        samtools view -h -@ 1 \
            --write-index \
            -o "relevant_reads.bam##idx##relevant_reads.bam.bai" \
            -M -L genotyping.sites.bed \
            ~{aligned_bam}
    >>>

    output {
        File relevant_reads     = "relevant_reads.bam"
        File relevant_reads_bai = "relevant_reads.bam.bai"
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:             4,
        mem_gb:                8,
        disk_gb:               375, # will use LOCAL SSD for speeding things up
        boot_disk_gb:          10,
        preemptible_tries:     0,
        max_retries:           1,
        docker:                "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                   select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:        select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:           select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:            select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                select_first([runtime_attr.docker, default_attr.docker])
    }
}

task ResetCLRBaseQual {
    input {
        File bam
        File bai

        Int arbitrary_bq
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 100 + 2*ceil(size(bam, "GB"))

    String prefix = "barbequed"

    command <<<
        set -eux

        python /usr/local/bin/reset_clr_bam_bq.py \
            -q ~{arbitrary_bq} \
            -p ~{prefix} \
            ~{bam}
        rm -f "~{prefix}.bai" "~{prefix}.bam.bai"
        samtools index "~{prefix}.bam"
    >>>

    output {
        File barbequed_bam = "~{prefix}.bam"
        File barbequed_bai = "~{prefix}.bam.bai"
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:             2,
        mem_gb:                8,
        disk_gb:               disk_size,
        boot_disk_gb:          10,
        preemptible_tries:     3,
        max_retries:           2,
        docker:                "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.34"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                   select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:        select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:           select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:            select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                select_first([runtime_attr.docker, default_attr.docker])
    }

}

task CheckCLRFingerprint {

    meta {
        description: "Uses Picard tool CheckFingerprint to verify if the samples in provided VCF and the CLR BAM arise from the same biological sample."
    }
    input {
        File aligned_bam
        File aligned_bai
        Int min_base_q = 0

        File fingerprint_vcf
        String vcf_sample_name

        File haplotype_map

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        haplotype_map: "table indicating reference sequence and auxillary file locations"
    }

    Int disk_size = 100 + ceil(size(aligned_bam, "GB"))
    String prefix = basename(aligned_bam, ".bam")

    command <<<
        set -eux

        java -jar /usr/picard/picard.jar \
            CheckFingerprint \
            INPUT=~{aligned_bam} \
            GENOTYPES=~{fingerprint_vcf} \
            EXPECTED_SAMPLE_ALIAS=~{vcf_sample_name} \
            HAPLOTYPE_MAP=~{haplotype_map} \
            OUTPUT=~{prefix} \
            MIN_BASE_QUAL=~{min_base_q}

        grep -v '^#' ~{prefix}.fingerprinting_summary_metrics | \
            grep -A1 READ_GROUP | \
            awk '
                {
                    for (i=1; i<=NF; i++)  {
                        a[NR,i] = $i
                    }
                }
                NF>p { p = NF }
                END {
                    for(j=1; j<=p; j++) {
                        str=a[1,j]
                        for(i=2; i<=NR; i++){
                            str=str" "a[i,j];
                        }
                        print str
                    }
                }' | \
            sed 's/ /\t/' \
            > metrics_map.txt

        mv ~{prefix}.fingerprinting_summary_metrics \
            ~{prefix}.fingerprinting_summary_metrics.txt
        mv ~{prefix}.fingerprinting_detail_metrics \
            ~{prefix}.fingerprinting_detail_metrics.txt
    >>>

    output {
        File summary_metrics = "~{prefix}.fingerprinting_summary_metrics.txt"
        File detail_metrics = "~{prefix}.fingerprinting_detail_metrics.txt"
        Map[String, String] metrics_map = read_map("metrics_map.txt")
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:             2,
        mem_gb:                4,
        disk_gb:               disk_size,
        boot_disk_gb:          10,
        preemptible_tries:     3,
        max_retries:           2,
        docker:                "us.gcr.io/broad-dsp-lrma/picard:lrfp-clr"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                   select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:        select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:           select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:            select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                select_first([runtime_attr.docker, default_attr.docker])
    }
}
