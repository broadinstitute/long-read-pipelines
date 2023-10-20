version 1.0

import "../../structs/Structs.wdl"

task Hiphase {

    meta {
        description: "Generates phased VCF. Note this runs fast so no need to parallize."
    }


    input {
        File bam
        File bai

        File unphased_snp_vcf
        File unphased_snp_tbi
        File unphased_sv_vcf
        File unphased_sv_tbi

        File ref_fasta
        File ref_fasta_fai
        String prefix

        Int memory = 200
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        RuntimeAttr? runtime_attr_override
    }

    Int bam_sz = ceil(size(bam, "GB"))
	Int disk_size = if bam_sz > 200 then 2*bam_sz else bam_sz + 200

    command <<<
        set -euxo pipefail

        # exclude records with no genotype information
        bcftools view -e 'GT="./."' ~{unphased_sv_vcf} -o ~{unphased_sv_vcf}.filtered.vcf


        hiphase \
        --threads 16 \
        --bam ~{bam} \
        --reference ~{ref_fasta} \
        --global-realignment-cputime 300 \
        --vcf ~{unphased_snp_vcf} \
        --output-vcf ~{unphased_snp_vcf}.phased.vcf.gz \
        --vcf ~{unphased_sv_vcf}.filtered.vcf \
        --output-vcf ~{unphased_sv_vcf}.phased.vcf.gz \
        --stats-file ~{prefix}.stats.csv \
        --blocks-file ~{prefix}.blocks.tsv \
        --summary-file ~{prefix}.summary.tsv

        tabix -p vcf ~{unphased_snp_vcf}.phased.vcf.gz
        tabix -p vcf ~{unphased_sv_vcf}.phased.vcf.gz
        
    >>>

    output {
        File phased_snp_vcf = "~{unphased_snp_vcf}.phased.vcf.gz"
        File phased_snp_vcf_tbi = "~{unphased_snp_vcf}.phased.vcf.gz.tbi"
        File phased_sv_vcf   = "~{unphased_sv_vcf}.phased.vcf.gz"
        File phased_sv_vcf_tbi = "~{unphased_sv_vcf}.phased.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          64,
        mem_gb:             200,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "hangsuunc/hiphase:0.7.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
