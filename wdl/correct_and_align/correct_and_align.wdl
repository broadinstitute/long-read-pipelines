workflow CorrectAndAlignWorkflow {
    String input_bam
    String output_prefix
    String base_image="kgarimella/pbtools@sha256:7bb26d8a142cd85999ed1a55a8b2b299dda6264709150eb968d90a694f6d5851"
    File ref_fasta

    call SplitSubreads {
        input:
            input_bam=input_bam,
            output_prefix=output_prefix,
            docker_image=base_image
    }

    scatter (subread_file in SplitSubreads.subread_files) {
        call Minimap2 {
            input:
                ref_fasta=ref_fasta,
                subread_file=subread_file,
                subread_aligned=basename(subread_file, ".bam") + ".aligned.bam",
                docker_image=base_image
        }

        call CCS {
            input:
                subread_file=subread_file,
                subread_ccs=basename(subread_file, ".bam") + ".ccs.bam",
                docker_image=base_image
        }

        call Minimap2CCS {
            input:
                ref_fasta=ref_fasta,
                subread_file=CCS.ccs,
                subread_aligned=basename(subread_file, ".bam") + ".ccs.aligned.bam",
                docker_image=base_image
        }
    }

    call MergeBams {
        input:
            bam_outs=Minimap2.aligned,
            merged_bam=basename(input_bam, ".bam") + ".aligned.merged.bam",
            docker_image=base_image
    }

    call MergeCCSBams {
        input:
            bam_outs=Minimap2CCS.aligned,
            merged_bam=basename(input_bam, ".bam") + ".ccs.aligned.merged.bam",
            docker_image=base_image
    }
}

task SplitSubreads {
    String input_bam
    String output_prefix
    String docker_image

    Int cpus = 1
    Int disk_size = 20
    Int num_reads_per_split = 5000

    command <<<
        java -jar /gatk.jar SplitSubreadsByZmw -I ${input_bam} -O ${output_prefix} -nr ${num_reads_per_split}
    >>>

    output {
        Array[File] subread_files = glob("*.bam")
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "2G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
    }
}

task Minimap2 {
    File ref_fasta
    File subread_file
    String subread_aligned
    String docker_image

    Int cpus = 1
    Int disk_size = 20

    command <<<
        pbmm2 align ${ref_fasta} ${subread_file} ${subread_aligned} --sort
    >>>

    output {
        File aligned = "${subread_aligned}"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "20G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
    }
}

task CCS {
    File subread_file
    String subread_ccs
    String docker_image

    Int cpus = 1
    Int disk_size = 20

    command <<<
        ccs --minLength 10000 --maxLength 16000 ${subread_file} ${subread_ccs}
    >>>

    output {
        File ccs = "${subread_ccs}"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "20G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
    }
}

task Minimap2CCS {
    File ref_fasta
    File subread_file
    String subread_aligned
    String docker_image

    Int cpus = 1
    Int disk_size = 20

    command <<<
        pbmm2 align ${ref_fasta} ${subread_file} ${subread_aligned} --sort --preset CCS
    >>>

    output {
        File aligned = "${subread_aligned}"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "20G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
    }
}

task MergeBams {
    Array[File] bam_outs
    String merged_bam
    String docker_image

    Int cpus = 1
    Int disk_size = 20

    command <<<
        set -euxo pipefail
        java -Xmx4g -jar /gatk.jar MergeSamFiles -I ${sep=" -I " bam_outs} -O ${merged_bam} -AS
    >>>

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "20G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
    }

    output {
        File merged = "${merged_bam}"
    }
}

task MergeCCSBams {
    Array[File] bam_outs
    String merged_bam
    String docker_image

    Int cpus = 1
    Int disk_size = 20

    command <<<
        set -euxo pipefail
        java -Xmx4g -jar /gatk.jar MergeSamFiles -I ${sep=" -I " bam_outs} -O ${merged_bam} -AS
    >>>

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "20G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
    }

    output {
        File merged = "${merged_bam}"
    }
}

