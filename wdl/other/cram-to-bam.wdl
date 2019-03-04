# THIS IS A STUB - NOT FULLY IMPLEMENTED OR TESTED

workflow CramToBamWorkflow {
    File ref_fasta
    File input_cram

    call CramToBam {
        input:
            ref_fasta=ref_fasta,
            input_bam=input_cram
    }
}

task CramToBam {
    File ref_fasta
    File input_cram
    Int disk_size = ceil(size(ref_fasta, "GB") + 2 * size(input_cram, "GB") + 10)

    command {
        set -xe

        samtools view -B ${input_cram}
    }

    output {
        File output_bam_bai="${input_bam}.bai"
        File output_bam_bai="${input_bam}.bai"
    }

    runtime {
        docker: "weisburd/base-image-for-str-tools@sha256:bc02c67c69bbd13165ef2cd83f7b2ec87814d99fc007070cc7cd8865b29998a9"
        disks: "local-disk ${disk_size} SSD"
    }
}
