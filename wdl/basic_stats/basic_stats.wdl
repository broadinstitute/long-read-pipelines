workflow BasicStatsWorkflow {
    String input_bam
    String output_prefix
    String base_image="kgarimella/pbeap-base@sha256:b2fdefb391c2d5537376c2cbf8d58700bab6f1e3e8ee14756a0087dabcc23009"

    call ReadLengths {
        input:
            input_bam=input_bam,
            output_prefix=output_prefix,
            docker_image=pbeap_base_image
    }
}

task ReadLengths {
    String input_bam
    String output_prefix
    String docker_image

    Int cpus = 1
    Int disk_size = 10

    command <<<
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        samtools view ${input_bam} | awk -F"[\t/]" '{ print $1, $2, $3, length($12) }' > ${output_prefix}.read_lengths.txt
    >>>

    output {
        File listing = "${output_prefix}.read_lengths.txt"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "1G"
        bootDiskSizeGb: 10
        disks: "local-disk ${disk_size} SSD"
    }
}
