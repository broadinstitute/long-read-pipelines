version 1.0

import "../../structs/Structs.wdl"


task LongphaseSNPs {

    meta {
        description: "Generates phased VCF. Note this runs fast so no need to parallize."
    }


    input {
        File bam
        File bai

        File unphased_snp_vcf
        File unphased_snp_tbi

        File ref_fasta
        File ref_fasta_fai

        String samplename
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        RuntimeAttr? runtime_attr_override
    }

    Int bam_sz = ceil(size(bam, "GB"))
	Int disk_size = if bam_sz > 200 then 2*bam_sz else bam_sz + 200

    command <<<
        set -euxo pipefail

        longphase phase \
        -s ~{unphased_snp_vcf} \
        -b ~{bam} \
        -r ~{ref_fasta} \
        -t 16 \
        -o ~{samplename} \
        --indels \
        --pb # or --pb for PacBio Hifi

        bcftools view ~{samplename}.vcf -Oz -o ~{samplename}.vcf.gz
        bcftools sort ~{samplename}.vcf.gz -O z -o ~{samplename}.sorted.vcf.gz
        tabix -p vcf ~{samplename}.sorted.vcf.gz
        
    >>>

    output {
        File phased_vcf = "~{samplename}.sorted.vcf.gz"
        File phased_vcf_tbi = "~{samplename}.sorted.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "hangsuunc/longphase:v1"
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

