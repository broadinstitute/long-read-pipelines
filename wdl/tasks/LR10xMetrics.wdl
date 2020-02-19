version 1.0

import "Structs.wdl"

workflow LR10xMetrics {
    input {
        File subreads
    }

    call Correct { input: subreads = subreads }

    output {
        File consensus = Correct.consensus
        File report = Correct.report
    }
}

task Correct {
    input {
        File subreads

        Int max_length = 21000
        Int min_passes = 2
        Int cpus = 4

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(subreads, "GB"))

    command <<<
        set -euxo pipefail

        ccs --max-length ~{max_length} --min-passes ~{min_passes} -j ~{cpus} ~{subreads} ccs_unmapped.bam
    >>>

    output {
        File consensus = "ccs_unmapped.bam"
        File report = "ccs_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             40,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.24"
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

