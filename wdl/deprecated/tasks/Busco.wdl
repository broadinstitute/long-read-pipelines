version 1.0

import "../../structs/Structs.wdl"

task Busco {
    input {
        File assembly
        String lineage # busco --list-datasets will list possible options

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        busco -m genome -i ~{assembly} -o busco -l ~{lineage}
        tar -czvf busco_output.tar.gz busco
    >>>

    output {
        File output_tar = "busco_output.tar.gz"
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:      8,
        mem_gb:         32,
        disk_gb:        50,
        boot_disk_gb:   25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "ezlabgva/busco:v4.1.1_cv1"
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