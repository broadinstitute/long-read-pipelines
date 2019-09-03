version 1.0

# TODO: describe purpose
task MergeBams {
    input {
        Array[File] aligned_shards
        String merged_name
        String docker
    }

    Int cpus = 2
    Int disk_size = 3*ceil(size(aligned_shards, "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx4g -jar /usr/local/bin/gatk.jar MergeSamFiles -I ~{sep=" -I " aligned_shards} -O ~{merged_name} -AS --CREATE_INDEX --USE_JDK_DEFLATER --USE_JDK_INFLATER --VALIDATION_STRINGENCY SILENT
    >>>

    output {
        File merged = "~{merged_name}"
        File merged_bai = basename(merged_name, ".bam") + ".bai"
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