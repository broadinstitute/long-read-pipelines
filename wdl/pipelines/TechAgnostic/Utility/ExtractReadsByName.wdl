version 1.0

import "../../../structs/Structs.wdl"

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/BAMutils.wdl" as BU


workflow ExtractReadsByName {
    meta {
        desciption: "Extract reads from a BAM file based on desired read names"
    }
    parameter_meta {
        qnames: "desired read names"
        prefix: "prefix of output files"
        mode: "We offer three modes of operations--0 for remote extraction at the expense of possible http-related failures, 1 for local extract at the expense of slower operations, 2 for trying remote first, and local if that fails"
    }
    input {
        File  bam
        File? bai
        File qnames
        String prefix

        Int mode
    }
    output {
        File selBam = output_bam
        File? selBai = output_bai
    }

    if (0>mode || 2<mode) {
        call Utils.StopWorkflow as InvalidInput { input: reason = "value for 'mode' must be one of [0, 1, 2]"}
    }

    if (0==mode) {
        call BU.ExtractReadsByNameRemote as remote {
            input: bam = bam, bai = bai, qnames = qnames, prefix = prefix
        }

        if (remote.failed) {
            call Utils.StopWorkflow { input: reason = "samtools streaming from bucket failed."}
        }
    }

    if (1==mode) {
        call BU.ExtractReadsByNameLocal as local {
            input: bam = bam, bai = bai, qnames = qnames, prefix = prefix
        }
    }

    if (2==mode) {
        call BU.ExtractReadsByNameRemote as FirstAttempt {
            input: bam = bam, bai = bai, qnames = qnames, prefix = prefix
        }

        if (FirstAttempt.failed) {
            call Utils.StopWorkflow as Unlucky { input: reason = "Samtools streaming from bucket failed."}
        }
        call BU.ExtractReadsByNameLocal as Insurance {
            input: bam = bam, bai = bai, qnames = qnames, prefix = prefix
        }
        File mode_2_bam = select_first([Insurance.selBam, FirstAttempt.selBam])
        if (defined(bai)) {
            File mode_2_bai = select_first([Insurance.selBai, FirstAttempt.selBai])
        }
    }

    File output_bam = select_first([remote.selBam, local.selBam, mode_2_bam])
    if (defined(bai)) {
        File output_bai = select_first([remote.selBai, local.selBai, mode_2_bai])
    }
}
