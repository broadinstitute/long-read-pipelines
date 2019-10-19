version 1.0

import "Structs.wdl"

# TODO: describe purpose
# Note: this task changes the incoming read group name
task CCS {
    input {
        File unmapped_shard
        String platform
        
        RuntimeAttr? runtime_attr_override
    }

    String ccs_shard_name = basename(unmapped_shard, ".bam") + ".ccs.corrected.bam"
    Int max_length = 21000
    Int min_passes = 2
    Int cpus = 4
    Int disk_size = 2*ceil(size(unmapped_shard, "GB"))

    command <<<
        set -euxo pipefail

        PLATFORM="~{platform}"

        if [ $PLATFORM == "ONT" ]
        then
            samtools view -H ~{unmapped_shard} | samtools view -b > ~{ccs_shard_name}

            echo "ZMWs input          (A)  : 0
ZMWs generating CCS (B)  : 0 (100.00%)
ZMWs filtered       (C)  : 0 (0.00%)

Exclusive ZMW counts for (C):
No usable subreads       : 0 (0.00%)
Below SNR threshold      : 0 (0.00%)
Lacking full passes      : 0 (0.00%)
Heteroduplexes           : 0 (0.00%)
Min coverage violation   : 0 (0.00%)
Draft generation error   : 0 (0.00%)
Draft above --max-length : 0 (0.00%)
Draft below --min-length : 0 (0.00%)
Lacking usable subreads  : 0 (0.00%)
CCS did not converge     : 0 (0.00%)
CCS below minimum RQ     : 0 (0.00%)
Unknown error            : 0 (0.00%)" > ccs_report.txt
        else
            ccs --max-length ~{max_length} --min-passes ~{min_passes} -j ~{cpus} ~{unmapped_shard} ~{ccs_shard_name}
        fi
    >>>

    output {
        File ccs_shard = "~{ccs_shard_name}"
        File ccs_report = "ccs_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          "~{cpus}", 
        mem_gb:             40, 
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.18"
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
