version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/PBUtils.wdl" as PB
import "../../../tasks/Utility/BAMutils.wdl" as BU

import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/QC/SampleLevelAlignedMetrics.wdl" as COV

workflow Work {
    meta {
        description: "Merge a sample's (potential) multiple SMRT-/flow-cells data, and collect alignment metrics."
    }
    parameter_meta {
        bams_suspected_to_contain_dup_record: "Some ONT output files from basecall dirs have a strange duplicate issue."
        bed_to_compute_coverage: "BED file holding regions-of-interest for computing coverage over."
        bed_descriptor: "Description of the BED file, will be used in the file name so be careful naming things"
    }
    input {
        String gcs_out_dir

        # sample specific
        String sample_name
        Array[File] aligned_bams
        Array[File] aligned_bais

        Boolean is_ont
        Boolean bams_suspected_to_contain_dup_record

        File? bed_to_compute_coverage
        String? bed_descriptor
    }
    output {
        File aligned_bam = FinalizeBam.gcs_path
        File aligned_bai = FinalizeBai.gcs_path
        File? aligned_pbi = FinalizePbi.gcs_path

        Float coverage = AlignmentMetrics.coverage
        File coverage_per_chr = FinalizeMosdepSumm.gcs_path
        File? bed_cov_summary = FinalizeRegionalCoverage.gcs_path

        Map[String, Float] alignment_metrics = AlignmentMetrics.reads_stats
    }
    if (defined(bed_to_compute_coverage)) {
        if (!defined(bed_descriptor)) {
            call Utils.StopWorkflow { input: reason = "Must provied descriptive name of the BED file if the file is provided."}
        }
    }

    String outdir = sub(gcs_out_dir, "/$", "") + "/alignments"

    ###########################################################
    # some input validation
    scatter (pair in zip(aligned_bams, aligned_bais)) {
        call Utils.InferSampleName {input: bam = pair.left, bai = pair.right}
    }
    call Utils.CheckOnSamplenames {input: sample_names = InferSampleName.sample_name}
    if (InferSampleName.sample_name[0] != sample_name) {
        call Utils.StopWorkflow as SM_mismatch { input: reason = "Provided sample name and those encoded in the BAM(s) don't match."}
    }

    ###########################################################
    # gather across (potential multiple) input BAMs
    if (length(aligned_bams) > 1) {
        call BU.MergeBamsWithSamtools as MergeAllReads { input: bams = aligned_bams, out_prefix = sample_name }
    }

    File bam = select_first([MergeAllReads.merged_bam, aligned_bams[0]])
    File bai = select_first([MergeAllReads.merged_bai, aligned_bais[0]])

    ###########################################################
    # ont specific: sometimes there are duplicate reads
    if (is_ont && bams_suspected_to_contain_dup_record) {
        call Utils.DeduplicateBam as RemoveONTDuplicates {
            input: aligned_bam = bam, aligned_bai = bai
        }
    }

    ###########################################################
    # save bam and index
    File use_this_bam = select_first([RemoveONTDuplicates.corrected_bam, bam])
    File use_this_bai = select_first([RemoveONTDuplicates.corrected_bai, bai])

    call FF.FinalizeToFile as FinalizeBam { input: outdir = outdir, file = use_this_bam, name = "~{sample_name}.bam" }
    call FF.FinalizeToFile as FinalizeBai { input: outdir = outdir, file = use_this_bai, name = "~{sample_name}.bam.bai" }

    ###########################################################
    # pacbio specific index
    if (! is_ont) {
        call PB.PBIndex as PBIndexSampleReads { input: bam = use_this_bam }
        call FF.FinalizeToFile as FinalizePbi { input: outdir = outdir, file = PBIndexSampleReads.pbi, name = "~{sample_name}.bam.pbi" }
    }
    ###########################################################
    call COV.SampleLevelAlignedMetrics as AlignmentMetrics {
        input:
            aligned_bam = use_this_bam,
            aligned_bai = use_this_bai,
            bed_to_compute_coverage = bed_to_compute_coverage,
            bed_descriptor = bed_descriptor
    }
    call FF.FinalizeToFile as FinalizeMosdepSumm { input: outdir = outdir, file = AlignmentMetrics.coverage_per_chr }

    if (defined(bed_to_compute_coverage)) {
        call FF.FinalizeToFile as FinalizeRegionalCoverage { input:
            outdir = outdir, file = select_first([AlignmentMetrics.bed_cov_summary]),
            name = "~{sample_name}.mosdepth_coverage.coverage_over_bed.~{bed_descriptor}.summary.txt"
        }
    }
}
