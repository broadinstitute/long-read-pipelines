version 1.0

import "Structs.wdl"

task Collapse {
    input {
        File vcf
        File tbi
        File ref_fasta

        String keep = "first" # {first,maxqual,common}

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([vcf, tbi, ref_fasta], "GB")) + 1

    command <<<
        set -euxo pipefail

        zcat ~{vcf} | grep -v 'SVTYPE=cnv' | bgzip > merged.no_cnvs.vcf.gz
        tabix -p vcf merged.no_cnvs.vcf.gz

        truvari collapse \
            -i merged.no_cnvs.vcf.gz \
            -c ~{prefix}.truvari_collapsed.vcf \
            -f ~{ref_fasta} \
            -T 4 \
            -k ~{keep} | \
            bcftools sort /dev/stdin -o ~{prefix}.truvari.vcf.gz -O z

        tabix ~{prefix}.truvari.vcf.gz
    >>>

    output {
        File collapsed_vcf = "~{prefix}.truvari.vcf.gz"
        File collapsed_tbi = "~{prefix}.truvari.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             24,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-truvari:3.5.0"
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