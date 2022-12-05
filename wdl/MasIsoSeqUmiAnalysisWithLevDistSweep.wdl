version 1.0

import "tasks/Structs.wdl"
import "tasks/MasIsoSeqUmiAnalysisTasks.wdl" as MAS_UMI
import "tasks/Utils.wdl" as Utils
import "tasks/Longbow.wdl" as LONGBOW
import "tasks/Finalize.wdl" as FF

workflow MasIsoSeqUmiAnalysisWithLevDistSweep {
    meta {
        description : "This workflow runs an analysis on the performance of the `umi cover` algorithm by replacing the UMIs in reads by a known adapter, and then applying that algorithm."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File pre_extracted_array_element_ccs_bam
        File pre_extracted_array_element_clr_bam

        File eq_class_tsv

        String prefix = "analyzed_umi_output"

        String eq_class_tag = "eq"
        String gene_tag = "XG"

        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/MasIsoSeqUmiAnalysisWithLevDistSweep"

        # Set up our parameter sweep:
        Array[Int] ccs_lev_dists = [0,0,1,2,3,4,5]
        Array[Int] clr_lev_dists = [0,1,2,3,4,5,6]

        # Min SW Score for excluding 3' adapters:
        Int min_sw_score = 35

        String sample_name
        String out_dir_suffix = ""
    }

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_001_WdlExecutionStartTimestamp { input: }

    # Merge Aligned CCS and Reclaimed reads together:
    call Utils.MergeBams as t_002_MergeCcsAndClrReads {
        input:
            bams = [pre_extracted_array_element_ccs_bam, pre_extracted_array_element_clr_bam],
            prefix = prefix + ".all_reads"
    }

    call MAS_UMI.RestoreSegCoordsAndOriginalNameToReads as t_003_RestoreSegCoordsAndOriginalNameToReads {
        input:
            bam_file = t_002_MergeCcsAndClrReads.merged_bam,
            bam_index = t_002_MergeCcsAndClrReads.merged_bai,
            prefix = prefix + ".all_reads.names_restored"
    }

    call MAS_UMI.CopyEqClassInfoToTag as t_004_CopyEqClassInfoToTag {
        input:
            bam_file = t_003_RestoreSegCoordsAndOriginalNameToReads.bam,
            bam_index = t_003_RestoreSegCoordsAndOriginalNameToReads.bai,
            eq_class_tsv = eq_class_tsv,
            eq_class_tag = eq_class_tag,
            gene_tag = gene_tag,
            prefix = prefix + ".all_reads.names_restored.eq_class_assigned",
    }

    call MAS_UMI.ExtractOptimial3pAdapterToUmiTag as t_005_ExtractOptimial3pAdapterToUmiTag {
        input:
            bam_file = t_004_CopyEqClassInfoToTag.bam,
            bam_index = t_004_CopyEqClassInfoToTag.bai,
            min_sw_score = min_sw_score,
            prefix = prefix + ".all_reads.names_restored.eq_class_assigned.3p_adapter_as_umi",
    }

    scatter (i1 in range(length(ccs_lev_dists))) {
        Int ccs_lev_dist = ccs_lev_dists[i1]
        Int clr_lev_dist = clr_lev_dists[i1]

        call MAS_UMI.UmiCoverForThreePrimeAnalysis as t_006_UmiCoverForThreePrimeAnalysis {
            input:
                bam = t_005_ExtractOptimial3pAdapterToUmiTag.bam,
                umi_tag = "ZV",
                final_umi_tag = "XX",
                back_alignment_score_tag = "ZB",
                umi_corrected_tag = "XZ",
                eq_class_tag = eq_class_tag,
                gene_tag = gene_tag,
                max_ccs_length_diff = 9999999,
                max_clr_length_diff = 9999999,
                max_ccs_gc_diff = 10,
                max_clr_gc_diff = 10,
                max_ccs_edit_dist = ccs_lev_dist,
                max_clr_edit_dist = clr_lev_dist,
                prefix = prefix + ".all_reads.names_restored.eq_class_assigned.3p_adapter_as_umi.umi_cover_corrected.ccs_lev_~{ccs_lev_dist}_clr_lev_~{clr_lev_dist}",
                runtime_attr_override = object {mem_gb: 128},
        }

        call MAS_UMI.SplitCcsAndClrReads as t_007_SplitCcsAndClrReads {
            input:
                bam_file = t_006_UmiCoverForThreePrimeAnalysis.umi_corrected_bam,
                bam_index = t_006_UmiCoverForThreePrimeAnalysis.umi_corrected_bam_index,
                prefix = prefix + ".names_restored.eq_class_assigned.3p_adapter_as_umi.umi_cover_corrected.ccs_lev_~{ccs_lev_dist}_clr_lev_~{clr_lev_dist}",
        }

        call MAS_UMI.CreateSimpleCountMatrixForUmiAnalysis as t_008_CreateSimpleCountMatrixForUmiAnalysisCCS {
            input:
                bam_file = t_007_SplitCcsAndClrReads.ccs_bam,
                bam_index = t_007_SplitCcsAndClrReads.ccs_bai,
                eq_class_tsv = eq_class_tsv,
                prefix = prefix + ".names_restored.eq_class_assigned.3p_adapter_as_umi.umi_cover_corrected.ccs_lev_~{ccs_lev_dist}_clr_lev_~{clr_lev_dist}.ccs.simple_counts",
        }
        call MAS_UMI.CreateSimpleCountMatrixForUmiAnalysis as t_009_CreateSimpleCountMatrixForUmiAnalysisCLR {
            input:
                bam_file = t_007_SplitCcsAndClrReads.clr_bam,
                bam_index = t_007_SplitCcsAndClrReads.clr_bai,
                eq_class_tsv = eq_class_tsv,
                prefix = prefix + ".names_restored.eq_class_assigned.3p_adapter_as_umi.umi_cover_corrected.ccs_lev_~{ccs_lev_dist}_clr_lev_~{clr_lev_dist}.clr.simple_counts",
        }
        call MAS_UMI.CreateSimpleCountMatrixForUmiAnalysis as t_010_CreateSimpleCountMatrixForUmiAnalysisAll {
            input:
                bam_file = t_006_UmiCoverForThreePrimeAnalysis.umi_corrected_bam,
                bam_index = t_006_UmiCoverForThreePrimeAnalysis.umi_corrected_bam_index,
                eq_class_tsv = eq_class_tsv,
                prefix = prefix + ".names_restored.eq_class_assigned.3p_adapter_as_umi.umi_cover_corrected.ccs_lev_~{ccs_lev_dist}_clr_lev_~{clr_lev_dist}.all.simple_counts",
        }
    }

    ################################################################################

    String outdir = sub(gcs_out_root_dir, "/$", "") + out_dir_suffix + "/" + sample_name + "/" + t_001_WdlExecutionStartTimestamp.timestamp_string
    File keyfile = select_first(t_008_CreateSimpleCountMatrixForUmiAnalysisCCS.simple_counts_tsv)

    call FF.FinalizeToDir as t_011_FinalizeNonSweptFiles {
            input:
                files = [
                    t_003_RestoreSegCoordsAndOriginalNameToReads.bam,
                    t_003_RestoreSegCoordsAndOriginalNameToReads.bai,

                    t_004_CopyEqClassInfoToTag.bam,
                    t_004_CopyEqClassInfoToTag.bai,

                    t_005_ExtractOptimial3pAdapterToUmiTag.bam,
                    t_005_ExtractOptimial3pAdapterToUmiTag.bai,
                    t_005_ExtractOptimial3pAdapterToUmiTag.ccs_levs_pickle,
                    t_005_ExtractOptimial3pAdapterToUmiTag.clr_levs_pickle,
                    t_005_ExtractOptimial3pAdapterToUmiTag.ccs_sws_pickle,
                    t_005_ExtractOptimial3pAdapterToUmiTag.clr_sws_pickle,
                    t_005_ExtractOptimial3pAdapterToUmiTag.rejected_bam_no_threep,
                    t_005_ExtractOptimial3pAdapterToUmiTag.rejected_bam_low_ssw_score,
                ],
                outdir = outdir + "/base_prepped_reads",
                keyfile = keyfile
     }

    scatter (i2 in range(length(ccs_lev_dists))) {
        Int ccsls = ccs_lev_dists[i2]
        Int clrld = clr_lev_dists[i2]

        call FF.FinalizeToDir as t_012_FinalizeSweptRefAndSt2Comparisons {
            input:
                files = [
                    t_006_UmiCoverForThreePrimeAnalysis.umi_corrected_bam[i2],
                    t_006_UmiCoverForThreePrimeAnalysis.umi_corrected_bam_index[i2],
                    t_006_UmiCoverForThreePrimeAnalysis.failed_umi_correction_bam[i2],
                    t_006_UmiCoverForThreePrimeAnalysis.cached_read_loci[i2],

                    t_007_SplitCcsAndClrReads.ccs_bam[i2],
                    t_007_SplitCcsAndClrReads.ccs_bai[i2],
                    t_007_SplitCcsAndClrReads.clr_bam[i2],
                    t_007_SplitCcsAndClrReads.clr_bai[i2],

                    t_008_CreateSimpleCountMatrixForUmiAnalysisCCS.simple_counts_tsv[i2],
                    t_009_CreateSimpleCountMatrixForUmiAnalysisCLR.simple_counts_tsv[i2],
                    t_010_CreateSimpleCountMatrixForUmiAnalysisAll.simple_counts_tsv[i2],
                ],
                outdir = outdir + "/ccs_lev_dist_~{ccsls}_clr_lev_dist_~{clrld}",
                keyfile = keyfile
        }
    }

    output {
        Array[File] all_reads_for_counting_bam = t_006_UmiCoverForThreePrimeAnalysis.umi_corrected_bam
        Array[File] all_reads_for_counting_bai = t_006_UmiCoverForThreePrimeAnalysis.umi_corrected_bam_index

        Array[File] ccs_reads_for_counting_bam = t_007_SplitCcsAndClrReads.ccs_bam
        Array[File] ccs_reads_for_counting_bai = t_007_SplitCcsAndClrReads.ccs_bai

        Array[File] clr_reads_for_counting_bam = t_007_SplitCcsAndClrReads.clr_bam
        Array[File] clr_reads_for_counting_bai = t_007_SplitCcsAndClrReads.clr_bai

        Array[File] ccs_simple_ccs_counts = t_008_CreateSimpleCountMatrixForUmiAnalysisCCS.simple_counts_tsv
        Array[File] clr_simple_ccs_counts = t_009_CreateSimpleCountMatrixForUmiAnalysisCLR.simple_counts_tsv
        Array[File] all_simple_ccs_counts = t_010_CreateSimpleCountMatrixForUmiAnalysisAll.simple_counts_tsv
    }
}
