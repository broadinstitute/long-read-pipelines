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

import "Utils.wdl" as Utils
import "PrepareData.wdl" as PD
import "ProcessReads.wdl" as PR
import "CallSVs.wdl" as CallSV
import "AssembleMT.wdl" as ASMT
import "AlignedMetrics.wdl" as AM
import "UnalignedMetrics.wdl" as UM
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
        File dbsnp
        File dbsnp_tbi
        File metrics_locus

        String gcs_output_dir

        Boolean? sample_is_female
    }

    String outdir = sub(sub(gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")

    scatter (gcs_dir in gcs_dirs) {
        call PD.PrepareData as PrepareData { input: gcs_dir = gcs_dir, sample_name = sample_name }

        String platform = PrepareData.run_info['PL']

        scatter (manifest_chunk in PrepareData.manifest_chunks) {
            call PR.ProcessReads as ProcessReads {
                input:
                    manifest_chunk = manifest_chunk,
                    ref_fasta = ref_fasta,
                    run_info = PrepareData.run_info
            }
        }

        call PR.FastqsToUnmappedBam as FastqsChunkToUnmappedBam {
            input:
                manifest_chunks = PrepareData.manifest_chunks,
                SM=PrepareData.run_info['SM'],
                ID=PrepareData.run_info['ID'],
                PL=PrepareData.run_info['PL']
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
                unaligned_bam = FastqsChunkToUnmappedBam.unmapped,

                per = "flowcell",
                type = "raw",
                label = PrepareData.run_info['ID'],

                gcs_output_dir = gcs_output_dir
        }

        call AM.AlignedMetrics as PerFlowcellAlignedMetrics {
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
                type = "raw",
                label = PrepareData.run_info['ID'],
                gcs_output_dir = gcs_output_dir
        }

        call ASMT.AssembleMT as AssembleMT {
            input:
                aligned_bam   = MergeAligned.merged,
                aligned_bai   = MergeAligned.merged_bai,

                platform      = PrepareData.run_info['PL'],

                ref_fasta     = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_dict      = ref_dict,
                mt_chr_name   = mt_chr_name,

                prefix        = PrepareData.run_info['PU'] + ".mt"
        }

        call FF.FinalizeToDir as FinalizeMT {
            input:
                files = [
                    AssembleMT.contigs_fasta,
                    AssembleMT.mt_aligned_bam,
                    AssembleMT.mt_aligned_bai,
                    AssembleMT.mt_calls
                ],
                outdir = outdir + "/metrics/per_flowcell/" + PrepareData.run_info['PU'] + "/mt/"
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

#    call PR.ValidateBam as ValidateAllRemaining {
#        input:
#            input_bam = MergeAllRemaining.merged,
#    }

    call CallSV.PBSV as PBSV {
        input:
            bam = MergeAllAligned.merged,
            bai = MergeAllAligned.merged_bai,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            tandem_repeat_bed = tandem_repeat_bed,
            prefix = basename(MergeAllAligned.merged, ".bam")
    }

    call CallSV.Sniffles as Sniffles {
        input:
            bam = MergeAllAligned.merged,
            bai = MergeAllAligned.merged_bai,
            prefix = basename(MergeAllAligned.merged, ".bam")
    }

    call CallSV.SVIM as SVIM {
        input:
            bam = MergeAllAligned.merged,
            bai = MergeAllAligned.merged_bai,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            prefix = basename(MergeAllAligned.merged, ".bam")
    }

    Array[String?] platform_gather = platform
    if ("PACBIO" == select_first(platform_gather)) {
        call PG.Peregrine as Peregrine {
            input:
                ref_fasta = ref_fasta,
                bam = MergeAllAligned.merged,
                sample_name = sample_name,
                output_prefix = basename(MergeAllAligned.merged, ".bam")
        }

        call DV.DeepVariant as DeepVariant {
            input:
                bam = MergeAllAligned.merged,
                bai = MergeAllAligned.merged_bai,
                ref_fasta = ref_fasta,
                ref_fai = ref_fasta_fai,
                model_class = "PACBIO",
                output_prefix = basename(MergeAllAligned.merged, ".bam")
        }

        call GATKBP.GATKBestPraciceForLR as GATKLR {
            input:
                input_bam = MergeAllAligned.merged,
                sample_is_female = sample_is_female,

                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_fai,
                ref_dict = ref_dict,

                output_prefix = basename(MergeAllAligned.merged, ".bam")
        }

        call FF.FinalizeToDir as FinalizePeregrineAssembly {
            input:
                files = [ Peregrine.final_fa, Peregrine.paf ],
                outdir = outdir + "/assembly"
        }

        call FF.FinalizeToDir as FinalizePeregrineCalls {
            input:
                files = [ Peregrine.variants ],
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
            files = [ PBSV.vcf, PBSV.tbi,
                      Sniffles.vcf, Sniffles.tbi,
                      SVIM.vcf, SVIM.tbi ],
            outdir = outdir + "/variants"
    }

    call FF.FinalizeToDir as FinalizeAlignedBams {
        input:
            files = [ MergeAllAligned.merged, MergeAllAligned.merged_bai ],
            outdir = outdir + "/alignments"
    }

    call FF.FinalizeToDir as FinalizeRemainingBams {
        input:
            files = [ MergeAllRemaining.merged, MergeAllRemaining.merged_bai ],
            outdir = outdir + "/alignments"
    }
}
