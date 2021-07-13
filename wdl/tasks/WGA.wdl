version 1.0

import "Structs.wdl"


task Minimap2 {
    input {
        File ref1
        File ref2

        String prefix

        String? map_preset = "asm5"
        RuntimeAttr? runtime_attr_override
    }

    output {
        File wga_bam = "~{prefix}.bam"
        File wga_bai = "~{prefix}.bam.bai"
    }

    command <<<
        minimap2 -ax ~{map_preset} -t ~{runtime_attr.cpu_cores} "~{ref1}" "~{ref2}" \
            | samtools sort -@ ~{runtime_attr.cpu_cores} -O bam -o "~{prefix}.bam"

        samtools index "~{prefix}.bam"
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-wga:0.1.0"
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


task Nucmer {
    input {
        File ref1
        File ref2

        String prefix
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        ref1: "Reference genome FASTA (has to be uncompressed)"
        ref2: "Query genome FASTA (has to be uncompressed)"

        prefix: "Output prefix. This task produces the original delta and a filtered file"
    }

    output {
        File delta = "~{prefix}.delta"
        File filtered_delta = "~{prefix}.filtered.delta"
    }

    command <<<
        set -euxo pipefail

        nucmer "~{ref1}" "~{ref2}" -p "~{prefix}"
        delta-filter -q "~{prefix}.delta" > "~{prefix}.filtered.delta"
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-wga:0.1.0"
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
