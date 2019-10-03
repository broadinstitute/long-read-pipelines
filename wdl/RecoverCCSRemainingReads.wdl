version 1.0

import "Structs.wdl"

# TODO: describe purpose
task RecoverCCSRemainingReads {
    input {
        File unmapped_shard
        File ccs_shard
        
        RuntimeAttr? runtime_attr_override
    }

    String remaining_shard_name = basename(unmapped_shard, ".bam") + ".ccs.remaining.bam"
    Int disk_size = 2*ceil(size(unmapped_shard, "GB"))

    command <<<
        set -euxo pipefail

        java -Dsamjdk.compression_level=0 -Xmx4g -jar /usr/local/bin/gatk.jar RecoverUncorrectedReads -I ~{unmapped_shard} -C ~{ccs_shard} -O ~{remaining_shard_name} -DF WellformedReadFilter --use-jdk-deflater --use-jdk-inflater
    >>>

    output {
        File remaining_shard = "~{remaining_shard_name}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2, 
        mem_gb:             20, 
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "kgarimella/lr-align:0.01.18"
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
