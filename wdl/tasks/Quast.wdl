version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.12/wdl/tasks/Structs.wdl"

task Quast {
    input {
        File ref
        Array[File] assemblies

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(ref, "GB") + size(assemblies, "GB"))

    command <<<
        set -euxo pipefail

        quast --no-icarus -r ~{ref}  ~{sep=' ' assemblies}

        tar czf quast_results.tgz quast_results/
    >>>

    output {
        File results = "quast_results.tgz"
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:      2,
        mem_gb:         4,
        disk_gb:        disk_size,
        boot_disk_gb:   10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-quast:0.1.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:        select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:     select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:     select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:        select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:         select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:             select_first([runtime_attr.docker, default_attr.docker])
    }
}