version 1.0

import "../../structs/Structs.wdl"

task HaploTagBam {
    meta {
        description: "Uses whatshap to haplotag a BAM."
    }

    parameter_meta {
        to_tag_bam: {description: "BAM to be haplotagged", localization_optional: true}
        phased_vcf: "VCF holding phased small variants to be used for the tagging"
    }

    input {
        File to_tag_bam
        File to_tag_bai
        File ref_fasta
        File ref_fasta_fai
        File phased_vcf
        File phased_tbi

        RuntimeAttr? runtime_attr_override
    }

    output {
        File tagged_bam = "~{prefix}.whatshap-haplotagged.bam"
        File tagged_bai = "~{prefix}.whatshap-haplotagged.bam.bai"
    }

    String prefix = basename(to_tag_bam, ".bam")
    String vcf_prefix = basename(phased_vcf, ".vcf.gz")

    Int disk_size = 10 + 2*ceil(size([to_tag_bam, phased_vcf], "GiB"))

    String local_bam = "/cromwell_root/~{prefix}.bam"
    String local_bai = "/cromwell_root/~{prefix}.bam.bai"

    command <<<
        set -eux

        time gcloud storage cp ~{to_tag_bam} ~{local_bam}
        mv ~{to_tag_bai} ~{local_bai}

        mv ~{phased_vcf} ~{vcf_prefix}.vcf.gz
        mv ~{phased_tbi} ~{vcf_prefix}.vcf.gz.tbi
        ls

        whatshap haplotag \
            -o ~{prefix}.whatshap-haplotagged.bam \
            --reference ~{ref_fasta} \
            --skip-missing-contigs \
            --output-threads=3 \
            ~{vcf_prefix}.vcf.gz \
            ~{local_bam}

        satmools index -@3 ~{prefix}.whatshap-haplotagged.bam
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-whatshap:2.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Phase {
    meta {
        description: "Uses whatshap phase (small variant) VCF."
    }

    parameter_meta {
        unphased_vcf: "VCF holding phased small variants to be phased."
        bam: {description: "BAM used for genrating the unphased VCF.", localization_optional: true}
    }

    input {
        File bam
        File bai
        File ref_fasta
        File ref_fasta_fai
        File unphased_vcf
        File unphased_tbi

        RuntimeAttr? runtime_attr_override
    }

    String bam_prefix = basename(bam, ".bam")
    String vcf_prefix = basename(unphased_vcf, ".vcf.gz")

    output {
        File phased_vcf = "~{vcf_prefix}.whatshap-phased.vcf.gz"
        File phased_tbi = "~{vcf_prefix}.whatshap-phased.vcf.gz.tbi"
    }

    Int disk_size = 10 + 2*ceil(size([bam, unphased_vcf], "GiB"))

    String local_bam = "/cromwell_root/~{bam_prefix}.bam"
    String local_bai = "~{local_bam}.bai"

    command <<<
        set -eux

        time gcloud storage cp ~{bam} ~{local_bam}
        mv ~{bai} ~{local_bai}
        mv ~{unphased_vcf} ~{vcf_prefix}.vcf.gz
        mv ~{unphased_tbi} ~{vcf_prefix}.vcf.gz.tbi
        ls

        whatshap phase \
            --indels \
            -o ~{vcf_prefix}.whatshap-phased.vcf.gz \
            --reference ~{ref_fasta} \
            ~{vcf_prefix}.vcf.gz \
            ~{local_bam}

        if [[ ! -f ~{vcf_prefix}.whatshap-phased.vcf.gz.tbi ]]; then tabix -p vcf ~{vcf_prefix}.whatshap-phased.vcf.gz; fi
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-whatshap:2.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Stats {
    meta {
        description: "Uses whatshap to haplotag a BAM."
    }

    parameter_meta {
        phased_vcf: "VCF holding phased small variants"
    }

    input {
        File phased_vcf
        File phased_tbi

        RuntimeAttr? runtime_attr_override
    }

    output {
        File stats_tsv = "~{vcf_prefix}.whatshap-stats.tsv"
        File stats_gtf = "~{vcf_prefix}.whatshap-stats.gtf"
    }

    String vcf_prefix = basename(phased_vcf, ".vcf.gz")

    Int disk_size = 10 + 2*ceil(size(phased_vcf, "GiB"))

    command <<<
        set -eux

        mv ~{phased_vcf} ~{vcf_prefix}.vcf.gz
        mv ~{phased_tbi} ~{vcf_prefix}.vcf.gz.tbi

        # for visualization
        whatshap stats \
            --gtf=~{vcf_prefix}.whatshap-stats.gtf \
            ~{vcf_prefix}.vcf.gz &

        # for use with MultiQC
        whatshap stats \
            --tsv=~{vcf_prefix}.whatshap-stats.tsv \
            ~{phased_vcf} &

        wait
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-whatshap:2.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
