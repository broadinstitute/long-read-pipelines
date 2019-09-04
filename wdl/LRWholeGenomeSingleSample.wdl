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

workflow LRWholeGenomeSingleSample {
    input {
        Array[String] gcs_dirs

        String? sample_name

        File ref_fasta
        File ref_fasta_fai
        File ref_dict
        String mt_chr_name
    }

    String docker_align = "kgarimella/lr-align:0.01.15"
    String docker_asm = "kgarimella/lr-asm:0.01.06"

    scatter (gcs_dir in gcs_dirs) {
        call Utils.DetectRunInfo as udri {
            input:
                gcs_dir = gcs_dir,
                sample_name = sample_name,
        }

        call Utils.PrepareRun as upr {
            input:
                files = udri.files,
        }

        call SLR.ShardLongReads as slrslr {
            input:
                unmapped_bam = upr.unmapped_bam,
        }

        scatter (unmapped_shard in slrslr.unmapped_shards) {
            call CR.CCS as crccs {
                input:
                    unmapped_shard = unmapped_shard,
                    platform = udri.run_info['PL'],
            }

            call AR.Minimap2 as AlignCCS {
                input:
                    shard = crccs.ccs_shard,
                    ref_fasta = ref_fasta,
                    SM = udri.run_info['SM'],
                    ID = udri.run_info['ID'] + ".corrected",
                    PL = udri.run_info['PL'],
                    reads_are_corrected = true,
            }

            call RCCSRR.RecoverCCSRemainingReads as rccsrr {
                input:
                    unmapped_shard = unmapped_shard,
                    ccs_shard = crccs.ccs_shard,
            }

            call AR.Minimap2 as AlignRemaining {
                input:
                    shard = rccsrr.remaining_shard,
                    ref_fasta = ref_fasta,
                    SM = udri.run_info['SM'],
                    ID = udri.run_info['ID'] + ".remaining",
                    PL = udri.run_info['PL'],
                    reads_are_corrected = false,
            }
        }

        call MB.MergeBams as MergeCorrected {
            input:
                aligned_shards = AlignCCS.aligned_shard,
                merged_name="corrected.bam",
        }

        call MB.MergeBams as MergeRemaining {
            input:
                aligned_shards = AlignRemaining.aligned_shard,
                merged_name="remaining.bam",
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

    call MB.MergeBams as MergeAllRemaining {
        input:
            aligned_shards = MergeRemaining.merged,
            merged_name="all.remaining.bam",
    }

    call VB.ValidateBam as ValidateAllRemaining {
        input:
            input_bam = MergeAllRemaining.merged,
    }
}
