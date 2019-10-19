version 1.0

# Copyright Broad Institute, 2019
#
# About:
#   This WDL pipeline processes long read data from a single sample (which may be split among multiple PacBio
#   SMRTCells or Oxford Nanopore flowcells).  We perform a variety of tasks (including CCS error correction,
#   alignment, flowcell merging, de novo assembly, SNP-to-SV variant discovery, variant filtration, methylation
#   calling, and automated QC. The results of pipeline are intended to be "analysis-ready".
#
#   This pipeline is capable of processing any combination of the following datasets:
#     - PacBio raw CCS data
#     - PacBio raw CLR data
#     - Oxford Nanopore raw data
#
#   Data may be presented as a PacBio run directory, Oxford Nanopore run directory, or just a collection of
#   BAM/fastq files.  Generally the data type and sample information are autodetected, but can also be manually
#   overridden in the input JSON file.
#
#
# Description of inputs:
#   Required:
#       Array[String] gcs_dirs          - The GCS directories wherein the data is stored.
#       String sample_name              - The sample name to use in place of the autodetected name.
#       File ref_fasta                  - The reference genome to which reads should be aligned.
#       File ref_fasta_fai              - The .fai index for the reference genome.
#       File ref_dict                   - The sequence dictionary for the reference genome.
#       String mt_chr_name              - The name of the contig in ref_fasta representing the mitochondrion genome.
#       File tandem_repeat_bed          - BED file representing tandem repeats in the reference.
#
#
# Licensing:
#   This script is released under the WDL source code license (BSD-3) (see LICENSE in
#   https://github.com/broadinstitute/wdl). Note however that the programs it calls may be subject to different
#   licenses. Users are responsible for checking that they are authorized to run all programs before running
#   this script.

import "AlignReads.wdl" as AR
import "CallSV.wdl" as CallSV
import "CorrectReads.wdl" as CR
import "MergeBams.wdl" as MB
import "RecoverCCSRemainingReads.wdl" as RCCSRR
import "ShardLongReads.wdl" as SLR
import "Utils.wdl" as Utils
import "ValidateBam.wdl" as VB
import "AssembleReads.wdl" as ASM
import "AssembleMT.wdl" as ASMT
import "LRMetrics.wdl" as MET
import "Peregrine.wdl" as PG
import "DeepVariantLR.wdl" as DV
import "GATKBestPractice.wdl" as GATKBP
import "Finalize.wdl" as FF

workflow LRWholeGenomeSingleSample {
    input {
        Array[String] gcs_dirs

        String sample_name

        File ref_fasta
        File ref_fasta_fai
        File ref_dict
        String mt_chr_name

        File tandem_repeat_bed
        File ref_flat

        String gcs_output_dir

        Boolean? sample_is_female
    }

    String outdir = sub(sub(gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")

    scatter (gcs_dir in gcs_dirs) {
        call Utils.DetectRunInfo as DetectRunInfo {
            input:
                gcs_dir = gcs_dir,
                sample_name = sample_name,
        }
        String platform = DetectRunInfo.run_info['PL']

        call Utils.PrepareRun as PrepareRun {
            input:
                files = DetectRunInfo.files,
        }

        call SLR.ShardLongReads as ShardLongReads {
            input:
                unmapped_bam = PrepareRun.unmapped_bam,
        }

        scatter (unmapped_shard in ShardLongReads.unmapped_shards) {
            call CR.CCS as CCS {
                input:
                    unmapped_shard = unmapped_shard,
                    platform = DetectRunInfo.run_info['PL'],
            }

            call AR.Minimap2 as AlignCCS {
                input:
                    shard = CCS.ccs_shard,
                    ref_fasta = ref_fasta,
                    SM = DetectRunInfo.run_info['SM'],
                    ID = DetectRunInfo.run_info['ID'] + ".corrected",
                    PL = DetectRunInfo.run_info['PL'],
                    reads_are_corrected = true,
            }

            call RCCSRR.RecoverCCSRemainingReads as RecoverCCSRemainingReads {
                input:
                    unmapped_shard = unmapped_shard,
                    ccs_shard = CCS.ccs_shard,
            }

            call AR.Minimap2 as AlignRemaining {
                input:
                    shard = RecoverCCSRemainingReads.remaining_shard,
                    ref_fasta = ref_fasta,
                    SM = DetectRunInfo.run_info['SM'],
                    ID = DetectRunInfo.run_info['ID'] + ".remaining",
                    PL = DetectRunInfo.run_info['PL'],
                    reads_are_corrected = false,
            }
        }

        call MB.MergeBams as MergeCorrected {
            input:
                aligned_shards = AlignCCS.aligned_shard,
                merged_name="~{sample_name}.corrected.bam",
        }

        call MB.MergeBams as MergeRemaining {
            input:
                aligned_shards = AlignRemaining.aligned_shard,
                merged_name="~{sample_name}.remaining.bam",
        }

        call ASMT.AssembleMT as AssembleMT {
            input:
                corrected_bam = MergeCorrected.merged,
                corrected_bai = MergeCorrected.merged_bai,
                remaining_bam = MergeRemaining.merged,
                remaining_bai = MergeRemaining.merged_bai,

                platform      = DetectRunInfo.run_info['PL'],

                ref_fasta     = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_dict      = ref_dict,
                mt_chr_name   = mt_chr_name,

                prefix        = DetectRunInfo.run_info['PU'] + ".mt"
        }

        call MET.LRMetrics as PerFlowcellMetrics {
            input:
                unaligned_bam = PrepareRun.unmapped_bam,
                aligned_bam = MergeCorrected.merged,
                aligned_bai = MergeCorrected.merged_bai,
                ref_dict = ref_dict,
                ref_flat = ref_flat,
        }

        # finalize per-flowcell metrics
        call FF.FinalizeToDir as FinalizeMT {
            input:
                files = [
                    AssembleMT.contigs_fasta,
                    AssembleMT.aligned_bam,
                    AssembleMT.aligned_bai,
                    AssembleMT.calls
                ],
                outdir = gcs_output_dir + "/metrics/per_flowcell/" + DetectRunInfo.run_info['PU'] + "/mt/"
        }

        #Array[File] coverage_full_dist      = MosDepth.full_dist
        call FF.FinalizeToDir as FinalizeCoverageFullGlobalDist {
            input:
                files = PerFlowcellMetrics.coverage_full_dist,
                outdir = gcs_output_dir + "/metrics/per_flowcell/" + DetectRunInfo.run_info['PU'] + "/coverage/"
        }

        #Array[File] coverage_global_dist    = MosDepth.global_dist
        call FF.FinalizeToDir as FinalizeCoverageGlobalDist {
            input:
                files = PerFlowcellMetrics.coverage_global_dist,
                outdir = gcs_output_dir + "/metrics/per_flowcell/" + DetectRunInfo.run_info['PU'] + "/coverage/"
        }

        #Array[File] coverage_region_dist    = MosDepth.region_dist
        call FF.FinalizeToDir as FinalizeCoverageRegionDist {
            input:
                files = PerFlowcellMetrics.coverage_region_dist,
                outdir = gcs_output_dir + "/metrics/per_flowcell/" + DetectRunInfo.run_info['PU'] + "/coverage/"
        }

        #Array[File] coverage_regions        = MosDepth.regions
        call FF.FinalizeToDir as FinalizeCoverageRegions {
            input:
                files = PerFlowcellMetrics.coverage_regions,
                outdir = gcs_output_dir + "/metrics/per_flowcell/" + DetectRunInfo.run_info['PU'] + "/coverage/"
        }

        #Array[File] coverage_regions_csi    = MosDepth.regions_csi
        call FF.FinalizeToDir as FinalizeCoverageSummary {
            input:
                files = PerFlowcellMetrics.coverage_regions_csi,
                outdir = gcs_output_dir + "/metrics/per_flowcell/" + DetectRunInfo.run_info['PU'] + "/coverage/"
        }

        #Array[File] coverage_quantized_dist = MosDepth.quantized_dist
        call FF.FinalizeToDir as FinalizeCoverageQuantizedDist {
            input:
                files = PerFlowcellMetrics.coverage_quantized_dist,
                outdir = gcs_output_dir + "/metrics/per_flowcell/" + DetectRunInfo.run_info['PU'] + "/coverage/"
        }

        #Array[File] coverage_quantized      = MosDepth.quantized
        call FF.FinalizeToDir as FinalizeCoverageQuantized {
            input:
                files = PerFlowcellMetrics.coverage_quantized,
                outdir = gcs_output_dir + "/metrics/per_flowcell/" + DetectRunInfo.run_info['PU'] + "/coverage/"
        }

        #Array[File] coverage_quantized_csi  = MosDepth.quantized_csi
        call FF.FinalizeToDir as FinalizeCoverageQuantizedCsi {
            input:
                files = PerFlowcellMetrics.coverage_quantized_csi,
                outdir = gcs_output_dir + "/metrics/per_flowcell/" + DetectRunInfo.run_info['PU'] + "/coverage/"
        }

        call FF.FinalizeToDir as FinalizeAlignedMetrics {
            input:
                files = [
                    PerFlowcellMetrics.aligned_flag_stats,
                    PerFlowcellMetrics.aligned_np_hist,
                    PerFlowcellMetrics.aligned_range_gap_hist,
                    PerFlowcellMetrics.aligned_zmw_hist,
                    PerFlowcellMetrics.aligned_prl_counts,
                    PerFlowcellMetrics.aligned_prl_hist,
                    PerFlowcellMetrics.aligned_prl_nx,
                    PerFlowcellMetrics.aligned_prl_yield_hist,
                    PerFlowcellMetrics.aligned_rl_counts,
                    PerFlowcellMetrics.aligned_rl_hist,
                    PerFlowcellMetrics.aligned_rl_nx,
                    PerFlowcellMetrics.aligned_rl_yield_hist,
                    PerFlowcellMetrics.rna_metrics
                ],
                outdir = gcs_output_dir + "/metrics/per_flowcell/" + DetectRunInfo.run_info['PU'] + "/aligned/"
        }

        call FF.FinalizeToDir as FinalizeUnalignedMetrics {
            input:
                files = [
                    PerFlowcellMetrics.unaligned_flag_stats,
                    PerFlowcellMetrics.unaligned_np_hist,
                    PerFlowcellMetrics.unaligned_range_gap_hist,
                    PerFlowcellMetrics.unaligned_zmw_hist,
                    PerFlowcellMetrics.unaligned_prl_counts,
                    PerFlowcellMetrics.unaligned_prl_hist,
                    PerFlowcellMetrics.unaligned_prl_nx,
                    PerFlowcellMetrics.unaligned_prl_yield_hist,
                    PerFlowcellMetrics.unaligned_rl_counts,
                    PerFlowcellMetrics.unaligned_rl_hist,
                    PerFlowcellMetrics.unaligned_rl_nx,
                    PerFlowcellMetrics.unaligned_rl_yield_hist
                ],
                outdir = gcs_output_dir + "/metrics/per_flowcell/" + DetectRunInfo.run_info['PU'] + "/unaligned/"
        }
    }

    call MB.MergeBams as MergeAllCorrected {
        input:
            aligned_shards = MergeCorrected.merged,
            merged_name="~{sample_name}.corrected.bam",
    }

    call VB.ValidateBam as ValidateAllCorrected {
        input:
            input_bam = MergeAllCorrected.merged,
    }

    call MB.MergeBams as MergeAllRemaining {
        input:
            aligned_shards = MergeRemaining.merged,
            merged_name="~{sample_name}.remaining.bam",
    }

    call VB.ValidateBam as ValidateAllRemaining {
        input:
            input_bam = MergeAllRemaining.merged,
    }

    call CallSV.PBSV as PBSV {
        input:
            bam = MergeAllCorrected.merged,
            bai = MergeAllCorrected.merged_bai,
            ref_fasta = ref_fasta,
            ref_fai = ref_fasta_fai,
            tandem_repeat_bed = tandem_repeat_bed,
            output_prefix = basename(MergeAllCorrected.merged, ".bam")
    }

    call CallSV.CompressAndIndex as CompressAndIndexPBSV {
        input:
            vcf = PBSV.variants
    }

    call CallSV.Sniffles as Sniffles {
        input:
            bam = MergeAllCorrected.merged,
            bai = MergeAllCorrected.merged_bai,
            sample_name = sample_name,
            output_prefix = basename(MergeAllCorrected.merged, ".bam")
    }

    call CallSV.CompressAndIndex as CompressAndIndexSniffles {
        input:
            vcf = Sniffles.variants
    }

    call CallSV.SVIM as SVIM {
        input:
            bam = MergeAllCorrected.merged,
            bai = MergeAllCorrected.merged_bai,
            ref_fasta = ref_fasta,
            ref_fai = ref_fasta_fai
    }

    Array[String?] platform_gather = platform
    if ("PACBIO" == select_first(platform_gather)) {
        call PG.Peregrine as Peregrine {
            input:
                ref_fasta = ref_fasta,
                bam = MergeAllCorrected.merged,
                sample_name = sample_name,
                output_prefix = basename(MergeAllCorrected.merged, ".bam")
        }

        call CallSV.CompressAndIndex as CompressAndIndexPeregrine {
            input:
                vcf = Peregrine.variants
        }

        call DV.DeepVariant as DeepVariant {
            input:
                bam = MergeAllCorrected.merged,
                bai = MergeAllCorrected.merged_bai,
                ref_fasta = ref_fasta,
                ref_fai = ref_fasta_fai,
                model_class = "PACBIO",
                output_prefix = basename(MergeAllCorrected.merged, ".bam")
        }

        call GATKBP.GATKBestPraciceForLR as GATKLR {
            input:
                input_bam = MergeAllCorrected.merged,
                sample_is_female = sample_is_female,

                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_fai,
                ref_dict = ref_dict,

                output_prefix = basename(MergeAllCorrected.merged, ".bam")
        }

        call FF.FinalizeToDir as FinalizePeregrineAssembly {
            input:
                files = [ Peregrine.final_fa, Peregrine.paf ],
                outdir = outdir + "/assembly"
        }

        call FF.FinalizeToDir as FinalizePeregrineCalls {
            input:
                files = [ CompressAndIndexPeregrine.variants, CompressAndIndexPeregrine.variants_tbi ],
                outdir = outdir + "/variants"
        }

        call FF.FinalizeToDir as FinalizeDV {
            input:
                files = [ DeepVariant.gvcf, DeepVariant.gvcf_tbi, DeepVariant.vcf, DeepVariant.vcf_tbi ],
                outdir = outdir + "/variants"
        }

        call FF.FinalizeToDir as FinalizeGATK {
            input:
                files = [ GATKLR.out_pp_vcf, GATKLR.out_pp_vcf_index, GATKLR.output_gvcf, GATKLR.output_gvcf_index ],
                outdir = outdir + "/variants"
        }
    }

    call FF.FinalizeToDir as FinalizeSVs {
        input:
            files = [ CompressAndIndexPBSV.variants, CompressAndIndexPBSV.variants_tbi,
                      CompressAndIndexSniffles.variants, CompressAndIndexSniffles.variants_tbi ],
            outdir = outdir + "/variants"
    }

    call FF.FinalizeToDir as FinalizeCorrectedBams {
        input:
            files = [ MergeAllCorrected.merged, MergeAllCorrected.merged_bai ],
            outdir = outdir + "/alignments"
    }

    call FF.FinalizeToDir as FinalizeUncorrectedBams {
        input:
            files = [ MergeAllRemaining.merged, MergeAllRemaining.merged_bai ],
            outdir = outdir + "/alignments"
    }
}
