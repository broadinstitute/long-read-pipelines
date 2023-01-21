version 1.0

import "Structs.wdl"

task Standardize {
    input {
        File vcf
        File ref_fai

        String caller
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([vcf, ref_fai], "GB")) + 1

    command <<<
        set -euxo pipefail

        svtk standardize \
            --include-reference-sites \
            --contigs ~{ref_fai} \
            --prefix ~{prefix} ~{vcf} - ~{caller} | \
            bcftools sort /dev/stdin -o ~{prefix}.truvari.std.vcf.gz -O z

        tabix ~{prefix}.truvari.std.vcf.gz
    >>>

    output {
        File standardized_vcf = "~{prefix}.truvari.std.vcf.gz"
        File standardized_tbi = "~{prefix}.truvari.std.vcf.gz.tbi"
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