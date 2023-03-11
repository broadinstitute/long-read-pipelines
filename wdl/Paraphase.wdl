version 1.0

import "tasks/Utils.wdl"

workflow Paraphase {
    meta {
        desciption: ""
    }
    input {
        File bam
        File bai
        String? sample_name
    }
    parameter_meta {

    }

    if (!defined(sample_name)) {
        call Utils.InferSampleName { input: bam = bam, bai = bai }
    }
    String smid = select_first([sample_name, InferSampleName.sample_name])
    call PacBioParaphase { input: bam = bam, bai = bai, sample_name = smid }

    output {
        File Paraphase_summary_json = PacBioParaphase.summary_json
        File Paraphase_res_bam = PacBioParaphase.realigned_tagged_bam
        File Paraphase_res_bai = PacBioParaphase.realigned_tagged_bai
        Array[File] Paraphase_vcfs = PacBioParaphase.vcfs
    }
}

task PacBioParaphase {
    meta {
        desciption: ""
    }
    input {
        File bam
        File bai
        String sample_name
    }
    parameter_meta {

    }

    command <<<
        set -eux

        paraphase \
            -b ~{bam} \
            -o output_directory \
            -v

        ls -r output_directory/
    >>>
    output {
        # Array[File] output_dir = glob("output_directory/*")

        File summary_json = "output_directory/~{sample_name}.json"
        File realigned_tagged_bam = "output_directory/~{sample_name}_realigned_tagged.bam"
        File realigned_tagged_bai = "output_directory/~{sample_name}_realigned_tagged.bam.bai"
        Array[File] vcfs = glob("output_directory/*.vcf")
    }

    runtime {
        cpu:    32
        memory: "192 GiB"
        disks:  "local-disk 375 LOCAL"
        docker: "us.gcr.io/broad-dsp-lrma/lr-paraphase:1.1.3-0"
    }
}
