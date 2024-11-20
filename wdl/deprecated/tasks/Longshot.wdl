version 1.0

import "../../structs/Structs.wdl"

# performs Longshot algo on one particular chromosome
task Longshot {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        File? sites_vcf
        File? sites_vcf_tbi

        Boolean phase = true

        String chr

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(select_all([bam, bai, ref_fasta, ref_fasta_fai, sites_vcf, sites_vcf_tbi]), "GB"))
    String prefix = basename(bam, ".bam")

    command <<<
        set -x

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)
        SM=$(samtools view -H ~{bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g')

        touch ~{prefix}.longshot.~{chr}.vcf
        longshot ~{true='-v' false='' defined(sites_vcf)} ~{select_first([sites_vcf, ""])} \
            -F \
            -r ~{chr} \
            -s $SM \
            --bam ~{bam} \
            --ref ~{ref_fasta} \
            ~{true='' false='--no_haps' phase} \
            --out ~{prefix}.longshot.~{chr}.vcf

        bgzip ~{prefix}.longshot.~{chr}.vcf
        tabix -p vcf ~{prefix}.longshot.~{chr}.vcf.gz
    >>>

    output {
        File vcf = "~{prefix}.longshot.~{chr}.vcf.gz"
        File vcf_tbi = "~{prefix}.longshot.~{chr}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             128,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longshot:0.4.1"
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
