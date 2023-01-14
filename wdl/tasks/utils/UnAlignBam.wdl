version 1.0

import "BAMutils.wdl" as BU

import "../Utils.wdl" as Utils

workflow UnAlignBam {
    meta {
        description: "Experiment that is now irrelevant given that samtools reset works better."
    }

    parameter_meta {

    }
    input {
        File bam
        File bai
        Boolean use_local_ssd
    }

    output {
        File uBAM = result_uBAM
    }

    if (use_local_ssd) {
        call Utils.ComputeAllowedLocalSSD as SSD { input: intended_gb = ceil(size(bam, "GB"))}
    }

    call BU.Drop2304Alignments { input: bam = bam, bai = bai, num_ssds = SSD.numb_of_local_ssd }
    call Utils.CountBamRecords as Original { input: bam = Drop2304Alignments.filtered_bam }
    call BU.SortBamByQName as QNSort { input: bam = Drop2304Alignments.filtered_bam, num_ssds = SSD.numb_of_local_ssd }





    # todo: clean up code in this function, and change to 1m reads

    call BU.SplitQNameSortedNo2304Bam as Chop {
        input:
            bam = QNSort.qnsort_bam, num_lines_per_chunk = 100000, num_ssds = SSD.numb_of_local_ssd
    }






    scatter (chunk in Chop.split_bams) {
        call BU.DropAlignmentInfo { input: bam = chunk }
    }
    Array[File] results = DropAlignmentInfo.res

    if (1 < length(results)) { # merge when necessary
        call BU.MergeBams as MergeUnaligned { input: bams = DropAlignmentInfo.res, prefix = basename(bam, ".bam") + ".unaligned", input_is_aligned = false }
    }
    File result_uBAM = select_first([MergeUnaligned.merged_bam, select_first(results)])

    call Utils.CountBamRecords as UnalignedCnt { input: bam = result_uBAM}

    if (Original.num_records != UnalignedCnt.num_records) {
        call Utils.StopWorkflow { input: reason = "Input primary & unmapped count " + Original.num_records + " doesn't match with unaligned version: " + UnalignedCnt.num_records }
    }
}
