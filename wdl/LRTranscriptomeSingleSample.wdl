version 1.0

# Copyright Broad Institute, 2019
#
# About:
#   This WDL pipeline processes long read RNA data from a single sample (which may be split among multiple
#   PacBio SMRTCells).  We perform a variety of tasks (including CCS error correction, alignment, and
#   flowcell merging).
#
#   This pipeline is currently capable of processing any combination of the following datasets:
#     - PacBio raw CCS IsoSeq data
#
#   Data may be presented as a PacBio run directory or just a collection of BAM/fastq files.  Generally the
#   data type and sample information are autodetected, but can also be manually overridden in the input JSON file.
#
#
# Description of inputs:
#   Required:
#       Array[String] gcs_dirs          - The GCS directories wherein the data is stored.
#       File ref_fasta                  - The reference genome to which reads should be aligned.
#       File ref_fasta_fai              - The .fai index for the reference genome.
#       File ref_dict                   - The sequence dictionary for the reference genome.
#
#   Optional:
#       String sample_name              - The sample name to use in place of the autodetected name.
#
#
# Licensing:
#   This script is released under the WDL source code license (BSD-3) (see LICENSE in
#   https://github.com/broadinstitute/wdl). Note however that the programs it calls may be subject to different
#   licenses. Users are responsible for checking that they are authorized to run all programs before running
#   this script.

import "AlignReads.wdl" as AR
import "CorrectReads.wdl" as CR
import "MergeBams.wdl" as MB
import "RecoverCCSRemainingReads.wdl" as RCCSRR
import "ShardLongReads.wdl" as SLR
import "Utils.wdl" as Utils
import "ValidateBam.wdl" as VB

workflow LRTranscriptomeSingleSample {
    input {
        Array[String] gcs_dirs

        String? sample_name

        File ref_fasta
        File ref_fasta_fai
        File ref_dict
    }

    scatter (gcs_dir in gcs_dirs) {
        call Utils.DetectRunInfo as DetectRunInfo {
            input:
                gcs_dir = gcs_dir,
                sample_name = sample_name,
        }

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
                    reads_are_rna = true,
            }
        }

        call MB.MergeBams as MergeCorrected {
            input:
                aligned_shards = AlignCCS.aligned_shard,
                merged_name="corrected.bam",
        }
    }

    call MB.MergeBams as MergeAllCorrected {
        input:
            aligned_shards = MergeCorrected.merged,
            merged_name="all.corrected.bam",
    }

    call VB.ValidateBam as ValidateAllCorrected {
        input:
            input_bam = MergeAllCorrected.merged,
    }
}
