version 1.0

import "../utils/BAMutils.wdl" as BU
import "../Utils.wdl"

workflow AlignmentFlagMetrics {
    input {
        File aligned_bam
        File aligned_bai
        Map[String, Int] names_and_decimal_flags = {'Secondary': 256,
                                                    'Supplementary': 2048,
                                                    'Unmapped': 4}
    }

    call Utils.ComputeAllowedLocalSSD { input: intended_gb = ceil(10 + size(aligned_bam, "GiB")) }
    call BU.CountAlignmentRecordsByFlag {
        input:
            aligned_bam = aligned_bam, aligned_bai = aligned_bai,
            names_and_decimal_flags = names_and_decimal_flags,
            num_local_ssds = ComputeAllowedLocalSSD.numb_of_local_ssd,
    }

    # Float all_cnt = All.count
    output {
        Int aln_record_count = CountAlignmentRecordsByFlag.total_cnt
        Map [String, Float] flat_metrics = CountAlignmentRecordsByFlag.flag_pcts
        Int summed_percentages = CountAlignmentRecordsByFlag.summed_percentages
    }
}

