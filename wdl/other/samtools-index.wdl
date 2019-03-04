
workflow SamtoolsIndexWorkflow {
    File input_bam_or_cram

    call SamtoolsIndex {
        input: input_bam_or_cram=input_bam_or_cram
    }
}

task SamtoolsIndex {
    File input_bam_or_cram
    Int disk_size = ceil(size(input_bam_or_cram, "GB") + 10)

    command {
        set -xe

        samtools index ${input_bam_or_cram}
    }

    output {
        File output_index="${input_bam_or_cram}.*ai"
    }

    runtime {
        docker: "weisburd/base-image-for-str-tools@sha256:bc02c67c69bbd13165ef2cd83f7b2ec87814d99fc007070cc7cd8865b29998a9"
        disks: "local-disk ${disk_size} HDD"
    }
}

