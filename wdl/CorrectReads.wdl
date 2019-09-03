version 1.0

# TODO: describe purpose
# Note: this task changes the incoming read group name
task CCS {
    input {
        File unmapped_shard
        String platform
        String docker
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

    runtime {
        cpu: "~{cpus}"
        memory: "40G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}