version 1.0

import "../../structs/Structs.wdl"


task HiphaseSNPs {

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

        hiphase \
        --threads 16 \
        --bam ~{bam} \
        --vcf ~{unphased_snp_vcf} \
        --output-vcf ~{samplename}_phased_snp.vcf.gz \
        --reference ~{ref_fasta} \

        bcftools sort ~{samplename}_phased_snp.vcf.gz -O z -o ~{samplename}_phased_snp.sorted.vcf.gz
        tabix -p vcf ~{samplename}_phased_snp.sorted.vcf.gz
        
    >>>

    output {
        File phased_vcf = "~{samplename}_phased_snp.sorted.vcf.gz"
        File phased_vcf_tbi = "~{samplename}_phased_snp.sorted.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "hangsuunc/hiphase:1.3.0"
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



task HiphaseSVs {

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
        String samplename

        Int memory = 200
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        RuntimeAttr? runtime_attr_override
    }

    Int bam_sz = ceil(size(bam, "GB"))
	Int disk_size = if bam_sz > 200 then 2*bam_sz else bam_sz + 200

    command <<<
        set -euxo pipefail

        touch ~{bai}

        hiphase \
        --threads 16 \
        --bam ~{bam} \
        --reference ~{ref_fasta} \
        --global-realignment-cputime 300 \
        --vcf ~{unphased_snp_vcf} \
        --output-vcf ~{samplename}_phased_snp.vcf.gz \
        --vcf ~{unphased_sv_vcf} \
        --output-vcf ~{samplename}_phased_sv.vcf.gz \
        --haplotag-file ~{samplename}_phased_sv_haplotag.tsv \
        --stats-file ~{samplename}.stats.csv \
        --blocks-file ~{samplename}.blocks.tsv \
        --summary-file ~{samplename}.summary.tsv \
        --verbose

        bcftools sort ~{samplename}_phased_snp.vcf.gz -O z -o ~{samplename}_phased_snp.sorted.vcf.gz
        tabix -p vcf ~{samplename}_phased_snp.sorted.vcf.gz

        bcftools sort ~{samplename}_phased_sv.vcf.gz -O z -o ~{samplename}_phased_sv.sorted.vcf.gz
        tabix -p vcf ~{samplename}_phased_sv.sorted.vcf.gz
        
    >>>

    output {
        File phased_snp_vcf = "~{samplename}_phased_snp.sorted.vcf.gz"
        File phased_snp_vcf_tbi = "~{samplename}_phased_snp.sorted.vcf.gz.tbi"
        File phased_sv_vcf   = "~{samplename}_phased_sv.sorted.vcf.gz"
        File phased_sv_vcf_tbi = "~{samplename}_phased_sv.sorted.vcf.gz.tbi"
        File haplotag_file = "~{samplename}_phased_sv_haplotag.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "hangsuunc/hiphase:1.3.0"
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


task HiphaseSVTrgt {

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
        File unphased_trgt_vcf
        File unphased_trgt_tbi

        File ref_fasta
        File ref_fasta_fai
        String samplename

        Int memory = 64
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        RuntimeAttr? runtime_attr_override
    }

    Int bam_sz = ceil(size(bam, "GB"))
	Int disk_size = if bam_sz > 200 then 2*bam_sz else bam_sz + 200

    command <<<
        set -euxo pipefail

        touch ~{bai}

        hiphase \
        --threads 16 \
        --bam ~{bam} \
        --reference ~{ref_fasta} \
        --global-realignment-cputime 300 \
        --vcf ~{unphased_snp_vcf} \
        --output-vcf ~{samplename}_phased_snp.vcf.gz \
        --vcf ~{unphased_sv_vcf} \
        --output-vcf ~{samplename}_phased_sv.vcf.gz \
        --vcf ~{unphased_trgt_vcf} \
        --output-vcf ~{samplename}_phased_trgt.vcf.gz \
        --haplotag-file ~{samplename}_phased_sv_haplotag.tsv \
        --stats-file ~{samplename}.stats.csv \
        --blocks-file ~{samplename}.blocks.tsv \
        --summary-file ~{samplename}.summary.tsv \
        

        bcftools sort ~{samplename}_phased_snp.vcf.gz -O z -o ~{samplename}_phased_snp.sorted.vcf.gz
        tabix -p vcf ~{samplename}_phased_snp.sorted.vcf.gz

        bcftools sort ~{samplename}_phased_sv.vcf.gz -O z -o ~{samplename}_phased_sv.sorted.vcf.gz
        tabix -p vcf ~{samplename}_phased_sv.sorted.vcf.gz
        
        bcftools sort ~{samplename}_phased_trgt.vcf.gz -O z -o ~{samplename}_phased_trgt.sorted.vcf.gz
        tabix -p vcf ~{samplename}_phased_trgt.sorted.vcf.gz
    >>>

    output {
        File phased_snp_vcf = "~{samplename}_phased_snp.sorted.vcf.gz"
        File phased_snp_vcf_tbi = "~{samplename}_phased_snp.sorted.vcf.gz.tbi"
        File phased_sv_vcf   = "~{samplename}_phased_sv.sorted.vcf.gz"
        File phased_sv_vcf_tbi = "~{samplename}_phased_sv.sorted.vcf.gz.tbi"
        File phased_trgt_vcf = "~{samplename}_phased_trgt.sorted.vcf.gz"
        File phased_trgt_vcf_tbi = "~{samplename}_phased_trgt.sorted.vcf.gz.tbi"
        File haplotag_file = "~{samplename}_phased_sv_haplotag.tsv"
        
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "hangsuunc/hiphase:1.3.0"
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
