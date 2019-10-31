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

import "Utils.wdl" as Utils
import "PrepareData.wdl" as PD
import "ProcessReads.wdl" as PR
import "AlignedMetrics.wdl" as AM
import "UnalignedMetrics.wdl" as UM
import "Finalize.wdl" as FF

workflow LRTranscriptomeSingleSample {
    input {
        Array[String] gcs_dirs

        String sample_name

        File ref_fasta
        File ref_fasta_fai
        File ref_dict
        File ref_flat

        File dbsnp
        File dbsnp_tbi
        File metrics_locus

        String gcs_output_dir
    }

    String outdir = sub(sub(gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")

    scatter (gcs_dir in gcs_dirs) {
        call PD.PrepareData as PrepareData { input: gcs_dir = gcs_dir, sample_name = sample_name, num_reads_per_split = 2000000 }
        String platform = PrepareData.run_info['PL']

        scatter (manifest_chunk in PrepareData.manifest_chunks) {
            call PR.ProcessReads as ProcessReads {
                input:
                    manifest_chunk = manifest_chunk,
                    ref_fasta = ref_fasta,
                    run_info = PrepareData.run_info,
                    reads_are_rna = true
            }
        }

        call PR.MergeBams as MergeRemaining {
            input:
                shards = ProcessReads.remaining_shard,
                merged_name="~{sample_name}.~{PrepareData.run_info['ID']}.unaligned.bam",
        }

        call PR.MergeBams as MergeAligned {
            input:
                shards = ProcessReads.aligned_shard,
                merged_name="~{sample_name}.~{PrepareData.run_info['ID']}.aligned.bam",
        }

        call UM.UnalignedMetrics as PerFlowcellUnalignedMetrics {
            input:
                unaligned_bam = PrepareData.files[0],

                per = "flowcell",
                type = "raw",
                label = PrepareData.run_info['ID'],

                gcs_output_dir = gcs_output_dir
        }

        call AM.AlignedMetrics as PerFlowcellConsensusMetrics {
            input:
                aligned_bam = MergeAligned.merged,
                aligned_bai = MergeAligned.merged_bai,
                ref_fasta = ref_fasta,
                ref_dict = ref_dict,
                ref_flat = ref_flat,
                dbsnp_vcf = dbsnp,
                dbsnp_tbi = dbsnp_tbi,
                metrics_locus = metrics_locus,
                per = "flowcell",
                type = "consensus",
                label = PrepareData.run_info['ID'],
                gcs_output_dir = gcs_output_dir
        }
    }

    call PR.MergeBams as MergeAllAligned {
        input:
            shards = MergeAligned.merged,
            merged_name="~{sample_name}.corrected.bam",
    }

    call PR.ValidateBam as ValidateAllAligned {
        input:
            input_bam = MergeAllAligned.merged,
    }

    call PR.MergeBams as MergeAllRemaining {
        input:
            shards = MergeRemaining.merged,
            merged_name="~{sample_name}.remaining.bam",
    }

    call PR.ValidateBam as ValidateAllRemaining {
        input:
            input_bam = MergeAllRemaining.merged,
    }

    call AM.AlignedMetrics as PerSampleConsensusMetrics {
        input:
            aligned_bam = MergeAllAligned.merged,
            aligned_bai = MergeAllAligned.merged_bai,
            ref_fasta = ref_fasta,
            ref_dict = ref_dict,
            ref_flat = ref_flat,
            dbsnp_vcf = dbsnp,
            dbsnp_tbi = dbsnp_tbi,
            metrics_locus = metrics_locus,
            per = "sample",
            type = "consensus",
            label = sample_name,
            gcs_output_dir = gcs_output_dir
    }

    call FF.FinalizeToDir as FinalizeAlignedBams {
        input:
            files = [ MergeAllAligned.merged, MergeAllAligned.merged_bai ],
            outdir = outdir + "/alignments"
    }
}
