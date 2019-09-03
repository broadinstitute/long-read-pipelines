version 1.0

# TODO: describe purpose
task ShardLongReads {
    input {
        File unmapped_bam
        String docker
        Int? num_reads_per_split
    }

    Int cpus = 1
    Int disk_size = 4*ceil(size(unmapped_bam, "GB"))
    Int nr = select_first([num_reads_per_split, 10000])

    command <<<
        set -euxo pipefail

        java -Dsamjdk.compression_level=0 -jar /usr/local/bin/gatk.jar ShardLongReads -I ~{unmapped_bam} -nr ~{nr} -O ./ -DF WellformedReadFilter --use-jdk-deflater --use-jdk-inflater
    >>>

    output {
        Array[File] unmapped_shards = glob("*.bam")
    }

    runtime {
        cpu: "~{cpus}"
        memory: "2G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}