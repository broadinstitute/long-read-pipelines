version 1.0

# TODO: describe purpose
task RecoverCCSRemainingReads {
    input {
        File unmapped_shard
        File ccs_shard
        String docker
    }

    String remaining_shard_name = basename(unmapped_shard, ".bam") + ".ccs.remaining.bam"
    Int cpus = 2
    Int disk_size = 2*ceil(size(unmapped_shard, "GB"))

    command <<<
        set -euxo pipefail

        java -Dsamjdk.compression_level=0 -Xmx4g -jar /usr/local/bin/gatk.jar RecoverUncorrectedReads -I ~{unmapped_shard} -C ~{ccs_shard} -O ~{remaining_shard_name} -DF WellformedReadFilter --use-jdk-deflater --use-jdk-inflater
    >>>

    output {
        File remaining_shard = "~{remaining_shard_name}"
    }

    runtime {
        cpu: "~{cpus}"
        memory: "20G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}