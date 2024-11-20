version 1.0

import "../../structs/Structs.wdl"

task ListGenotypedVCFs {

    meta {
        description: "List all VCFs that have been genotyped"
    }
    parameter_meta {
        fingerprint_store: "A bucket that holds all fingerprinting VCFs"
    }

    input {
        String fingerprint_store
    }

    String bucket_dir = sub(fingerprint_store, "/$", "")

    command <<<
        set -eux

        gsutil ls -r ~{bucket_dir}/**.vcf.gz > all.vcfs.txt
    >>>

    output {
        File vcf_gs_paths = "all.vcfs.txt"
    }

    ###################
    runtime {
        cpu: 2
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 25
        docker:"us.gcr.io/broad-dsp-lrma/lr-basic:latest"
    }
}

task PickGenotypeVCF {

    meta {
        description: "Pick a subset of VCFs that have been genotyped"
    }
    parameter_meta {
        fingerprinting_vcf_gs_paths:    "A file holding GS paths to fingerprinting GT'ed VCFs"
        vcf_name: "an expression used for picking up VCFs, the filter will be applied to VCF names, a match will lead to the VCF to be included"
    }

    input {
        File fingerprinting_vcf_gs_paths
        String? vcf_name
    }

    Boolean filter = defined(vcf_name)

    command <<<
        set -eux

        if ~{filter}; then
            grep "~{vcf_name}$" ~{fingerprinting_vcf_gs_paths} > vcfs.txt
        else
            cp ~{fingerprinting_vcf_gs_paths} vcfs.txt
        fi
    >>>

    output {
        Array[String] vcfs = if filter then [read_string("vcfs.txt")] else read_lines("vcfs.txt")
    }

    ###################
    runtime {
        cpu: 2
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 25
        docker:"gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task FilterGenotypesVCF {

    meta {
        description: "Filter out unwanted genotypes from a VCF"
    }
    parameter_meta {
        fingerprint_vcf: "A VCF file that has been genotyped"
        filters: "An array of chromosome names to filter out when verifying fingerprints"
    }

    input {
        File fingerprint_vcf
        Array[String] filters = ['_random\\t', '_decoy\\t', '_alt\\t', '^chrUn', '^HLA', '^EBV']
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
        File ready_to_use_vcf = "fingerprint.fixed.vcf"
    }

    ###################
    runtime {
        cpu: 2
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 25
        preemptible:     3
        maxRetries:      2
        docker:"gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task ExtractGenotypingSites {

    meta {
        description: "Extract genotyping sites from a VCF"
    }
    parameter_meta {
        fingerprint_vcf: "A VCF file that has been genotyped"
    }

    input {
        File fingerprint_vcf
    }

    command <<<

        set -eux

        GREPCMD="grep"
        if [[ ~{fingerprint_vcf} =~ \.gz$ ]]; then
            GREPCMD="zgrep"
        fi

        "${GREPCMD}" -v "^#" ~{fingerprint_vcf} | \
            awk 'BEGIN {OFS="\t"} {print $1, $2-1, $2, $3}' \
            > genotyping.sites.bed
    >>>

    output {
        File sites = "genotyping.sites.bed"
    }

    ###################
    runtime {
        cpu: 2
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 25
        preemptible:     3
        maxRetries:      2
        docker:"gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task MergeGenotypingSites {

    meta {
        description: "Merge genotyping sites from multiple VCFs"
    }
    parameter_meta {
        all_sites: "An array of genotyping sites"
    }

    input {
        Array[File] all_sites
    }

    command <<<

        set -eux
        cat ~{sep=' ' all_sites} | sort | uniq > "genotyping.sites.union.bed"
    >>>

    output {
        File merged_sites = "genotyping.sites.union.bed"
    }

    ###################
    runtime {
        cpu: 2
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 25
        preemptible:     3
        maxRetries:      2
        docker:"gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task ExtractRelevantGenotypingReads {

    meta {
        description: "Based on genotyping (SNP) sites, extract reads that overlap those places"
    }
    parameter_meta {
        aligned_bam:{
            localization_optional: true
        }
    }

    input {
        File aligned_bam
        File aligned_bai

        File genotyping_sites_bed
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 50 + 2*ceil(size(aligned_bam, "GB")) + 2*ceil(size(genotyping_sites_bed, "GB")) + 2*ceil(size(aligned_bai, "GB"))

    command <<<

        set -eux

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

        samtools view -h -@ 2 \
            --write-index \
            -o "relevant_reads.bam##idx##relevant_reads.bam.bai" \
            -M -L ~{genotyping_sites_bed} \
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
        disk_gb:               disk_size,
        boot_disk_gb:          25,
        preemptible_tries:     0,
        max_retries:           1,
        docker:                "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                   select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"  # if SSD is too slow, revert to LOCAL
        bootDiskSizeGb:        select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:           select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:            select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                select_first([runtime_attr.docker, default_attr.docker])
    }
}

task ResetCLRBaseQual {

    meta {
        description: "Reset base quality scores to a fixed value"
    }
    parameter_meta {
        bam: "A BAM file"
        bai: "A BAI file"
        arbitrary_bq: "A fixed base quality score"
    }

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
        boot_disk_gb:          25,
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

task CheckFingerprint {

    meta {
        description: "Uses Picard tool CheckFingerprint to verify if the samples in provided VCF and BAM arise from the same biological sample"
    }
    parameter_meta {
        aligned_bam:{
            description:  "GCS path to aligned BAM file, supposed to be of the same sample as from the fingerprinting VCF",
            localization_optional: true
        }

        fingerprint_vcf:    "Fingerprint VCF file from local database; note that sample name must be the same as in BAM"
        vcf_sample_name:    "Sample name in VCF, possibly different from that in the BAM."
        haplotype_map:      "Happlotype map file for the reference build used. See https://bit.ly/3QyZbwt"
    }

    input {
        File aligned_bam
        File aligned_bai

        File fingerprint_vcf
        String vcf_sample_name

        File haplotype_map

        RuntimeAttr? runtime_attr_override
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
        boot_disk_gb:          25,
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

task CheckCLRFingerprint {

    meta {
        description: "Uses Picard tool CheckFingerprint to verify if the samples in provided VCF and the CLR BAM arise from the same biological sample."
    }
    parameter_meta {
        vcf_sample_name:    "Sample name in VCF, possibly different from that in the BAM."
        haplotype_map:      "Haplotype map file for the reference build used. See https://bit.ly/3QyZbwt"
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
        boot_disk_gb:          25,
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

task ReheaderFullGRCh38VCFtoNoAlt {

    meta {
        desciption:
        "Reheader the fingperint VCF that's generated with full GRCh38 reference to the no_alt header; project specific."
    }
    parameter_meta {
        full_GRCh38_vcf:    "Full GRCh38 VCF file."
    }

    input {
        File full_GRCh38_vcf
    }

    command <<<
        set -eux

        GREPCMD="grep"
        if [[ ~{full_GRCh38_vcf} =~ \.gz$ ]]; then
            GREPCMD="zgrep"
        fi
        "${GREPCMD}" -vF "_decoy,length=" ~{full_GRCh38_vcf} | \
            grep -vF "_alt,length=" | \
            grep -v "^##contig=<ID=HLA-" \
            > "reheadered.fp.vcf"
    >>>

    output {
        File reheadered_vcf = "reheadered.fp.vcf"
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
