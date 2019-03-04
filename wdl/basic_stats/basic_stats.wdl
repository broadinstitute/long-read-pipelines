# see https://github.com/Oshlack/STRetch/wiki/Running-STRetch

workflow BasicStatsWorkflow {
    String input_bam
    String output_prefix
    String pbeap_base_image="kgarimella/pbeap-base"

    call FileSize {
        input:
            input_bam=input_bam,
            output_prefix=output_prefix,
            docker_image=pbeap_base_image
    }
}

task FileSize {
    String input_bam
    String output_prefix
    String docker_image

    Int cpus = 1
    Int disk_size = 1

    command {
	ls -lah ${input_bam} > ${output_prefix}.listing.txt
    }

    output {
	File listing = "*.listing.txt"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "1G"
        disks: "local-disk ${disk_size} SSD"
    }
}
