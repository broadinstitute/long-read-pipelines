version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF

import "tasks/SampleLevelAlignedMetrics.wdl" as COV

workflow ONTFlowcellFromMultipleBasecalls {
    input {
        Array[File] aligned_bams
        Array[File] aligned_bais
        Boolean bams_suspected_to_contain_dup_record = true

        String flowcell

        File? bed_to_compute_coverage

        File ref_map_file

        String gcs_out_root_dir
    }

    parameter_meta {
        bams_suspected_to_contain_dup_record: "if the multiple basecalls are suspected to have duplicates amongst them"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTFlowcell/~{flowcell}"

    ## Merge & deduplicate
    if (length(aligned_bams) > 1) {
        call Utils.MergeBams as MergeAllReads { input: bams = aligned_bams, prefix = flowcell }
    }

    File bam = select_first([MergeAllReads.merged_bam, aligned_bams[0]])
    File bai = select_first([MergeAllReads.merged_bai, aligned_bais[0]])
    if (bams_suspected_to_contain_dup_record) {
        call Utils.DeduplicateBam as RemoveDuplicates {
            input: aligned_bam = bam, aligned_bai = bai, same_name_as_input = true
        }
    }
    File usable_bam = select_first([RemoveDuplicates.corrected_bam, bam])
    File usable_bai = select_first([RemoveDuplicates.corrected_bai, bai])

    # collect metrics
    call COV.SampleLevelAlignedMetrics as coverage {
        input:
            aligned_bam = usable_bam,
            aligned_bai = usable_bai,
            ref_fasta   = ref_map['fasta'],
            bed_to_compute_coverage = bed_to_compute_coverage
    }

    # Finalize data
    String dir = outdir + "/alignments"

    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = dir, file = usable_bam }
    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = dir, file = usable_bai }
    if (defined(bed_to_compute_coverage)) { call FF.FinalizeToFile as FinalizeRegionalCoverage { input: outdir = dir, file = select_first([coverage.bed_cov_summary]) } }

    output {
        # Aligned BAM file
        File aligned_bam = FinalizeAlignedBam.gcs_path
        File aligned_bai = FinalizeAlignedBai.gcs_path

        # Unaligned read stats, please check entries in ONTBasecall

        # Aligned read stats
        File? bed_cov_summary = FinalizeRegionalCoverage.gcs_path

        Float aligned_num_reads = coverage.aligned_num_reads
        Float aligned_num_bases = coverage.aligned_num_bases
        Float aligned_frac_bases = coverage.aligned_frac_bases
        Float aligned_est_fold_cov = coverage.aligned_est_fold_cov

        Float aligned_read_length_mean = coverage.aligned_read_length_mean
        Float aligned_read_length_median = coverage.aligned_read_length_median
        Float aligned_read_length_stdev = coverage.aligned_read_length_stdev
        Float aligned_read_length_N50 = coverage.aligned_read_length_N50

        Float average_identity = coverage.average_identity
        Float median_identity = coverage.median_identity
    }
}