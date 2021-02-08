version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/AnnotateAdapters.wdl" as AA
import "tasks/Figures.wdl" as FIG
import "tasks/Finalize.wdl" as FF

workflow PBCCS10x {
    input {
        Array[File] bams
        File ref_map_file

        String participant_name
        File barcode_file
        Int num_shards = 300

        String? gcs_out_root_dir
    }

    parameter_meta {
        bams:             "GCS path to raw subreads or CCS data"
        ref_map_file:     "table indicating reference sequence and auxillary file locations"

        participant_name: "name of the participant from whom these samples were obtained"
        barcode_file:     "GCS path to the fasta file that specifies the expected set of multiplexing barcodes"
        num_shards:       "[default-valued] number of sharded BAMs to create (tune for performance)"

        gcs_out_root_dir: "[optional] GCS bucket to store the corrected/uncorrected reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    call Utils.GetDefaultDir { input: workflow_name = "PBCCS10x" }
    String outdir = sub(select_first([gcs_out_root_dir, GetDefaultDir.path]), "/$", "") + "/" + participant_name

    # scatter over all sample BAMs
    scatter (bam in bams) {
        File pbi = sub(bam, ".bam$", ".bam.pbi")

        call PB.GetRunInfo { input: bam = bam }
        String ID = GetRunInfo.run_info["PU"]

        # break one raw BAM into fixed number of shards
        call PB.ShardLongReads { input: unaligned_bam = bam, unaligned_pbi = pbi, num_shards = num_shards }

        # then perform correction on each of the shard
        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads }

            call AA.AnnotateAdapters { input: bam = CCS.consensus }

            call PB.Align as AlignCorrected {
                input:
                    bam         = AnnotateAdapters.annotated_bam,
                    ref_fasta   = ref_map['fasta'],
                    sample_name = participant_name,
                    map_preset  = "ISOSEQ"
            }

            call Utils.CountBamRecords as CountSubreadsInShard { input: bam = subreads }
            call Utils.CountBamRecords as CountCorrectedReadsInShard { input: bam = CCS.consensus }
            call Utils.CountBamRecords as CountAnnotatedReadsInShard { input: bam = AnnotateAdapters.annotated_bam }
        }

#        call C3.Cat as CountNumPassesInRun { input: files = CountNumPasses.num_passes, out = "num_passes.txt" }

        call Utils.Sum as CountSubreadsInFlowcell { input: ints = CountSubreadsInShard.num_records }
        call Utils.Sum as CountAnnotatedReadsInFlowcell { input: ints = CountAnnotatedReadsInShard.num_records }
        call Utils.Sum as CountCorrectedReadsInFlowcell { input: ints = CountCorrectedReadsInShard.num_records }

        call Utils.MergeBams as MergeAnnotated { input: bams = AnnotateAdapters.annotated_bam }
        call Utils.MergeBams as MergeCorrected { input: bams = AlignCorrected.aligned_bam }

        # compute alignment metrics
        call AM.AlignedMetrics as PerFlowcellMetrics {
            input:
                aligned_bam    = MergeCorrected.merged_bam,
                aligned_bai    = MergeCorrected.merged_bai,
                ref_fasta      = ref_map['fasta'],
                ref_dict       = ref_map['dict'],
                gcs_output_dir = outdir + "/metrics/per_flowcell/" + ID
        }

        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report, prefix = "~{participant_name}.~{ID}" }

        call FF.FinalizeToDir as FinalizeCCSReport {
            input:
                files = [ MergeCCSReports.report ],
                outdir = outdir + "/metrics/per_flowcell/" + ID + "/ccs_metrics"
        }
    }

#    call C3.Cat as CountNumPassesAll { input: files = CountNumPassesInRun.merged, out = "num_passes.txt" }

    call Utils.Sum as CountSubreads { input: ints = CountSubreadsInFlowcell.sum, prefix = "num_subreads" }
    call Utils.Sum as CountCorrectedReads { input: ints = CountCorrectedReadsInFlowcell.sum, prefix = "num_consensus" }
    call Utils.Sum as CountAnnotatedReads { input: ints = CountAnnotatedReadsInFlowcell.sum, prefix = "num_annotated" }

    # gather across (potential multiple) input raw BAMs
    if (length(bams) > 1) {
        call Utils.MergeBams as MergeAllCorrected { input: bams = MergeCorrected.merged_bam, prefix = "~{participant_name}.corrected" }
        call Utils.MergeBams as MergeAllAnnotated { input: bams = MergeAnnotated.merged_bam, prefix = "~{participant_name}.annotated" }

        call PB.MergeCCSReports as MergeAllCCSReports { input: reports = MergeCCSReports.report }
    }

    File ccs_bam = select_first([ MergeAllCorrected.merged_bam, MergeCorrected.merged_bam[0] ])
    File ccs_bai = select_first([ MergeAllCorrected.merged_bai, MergeCorrected.merged_bai[0] ])
    File ccs_report = select_first([ MergeAllCCSReports.report, MergeCCSReports.report[0] ])

    File annotated_bam = select_first([ MergeAllAnnotated.merged_bam, MergeAnnotated.merged_bam[0] ])
    File annotated_bai = select_first([ MergeAllAnnotated.merged_bai, MergeAnnotated.merged_bai[0] ])

    call Utils.GrepCountBamRecords as GrepAnnotatedReadsWithCBC {
        input:
            bam = annotated_bam,
            prefix = "num_annotated_with_cbc",
            regex = "CB:Z:[ACGT]"
    }

    call Utils.GrepCountBamRecords as GrepAnnotatedReadsWithCBCAndUniqueAlignment {
        input:
            bam = annotated_bam,
            samfilter = "-F 0x100",
            prefix = "num_annotated_with_cbc_and_unique_alignment",
            regex = "CB:Z:[ACGT]"
    }

    call Utils.BamToTable { input: bam = annotated_bam, prefix = "reads_aligned_annotated.table" }

    # compute alignment metrics
    call AM.AlignedMetrics as PerSampleMetrics {
        input:
            aligned_bam    = annotated_bam,
            aligned_bai    = annotated_bai,
            ref_fasta      = ref_map['fasta'],
            ref_dict       = ref_map['dict'],
            gcs_output_dir = outdir + "/metrics/combined/" + participant_name
    }

    ##########
    # Finalize
    ##########

    call FF.FinalizeToDir as FinalizeCorrectedReadCounts {
        input:
            files = [ CountSubreads.sum_file, CountAnnotatedReads.sum_file, CountCorrectedReads.sum_file,
                      GrepAnnotatedReadsWithCBC.num_records_file, GrepAnnotatedReadsWithCBCAndUniqueAlignment.num_records_file ],
            outdir = outdir + "/metrics/read_counts"
    }

#    call FF.FinalizeToDir as FinalizeNumPasses {
#        input:
#            files = [ CountNumPassesAll.merged ],
#            outdir = outdir + "/metrics/num_passes"
#    }

    call FF.FinalizeToDir as FinalizeBamTable {
        input:
            files = [ BamToTable.table ],
            outdir = outdir + "/metrics/bam_tables"
    }

    call FF.FinalizeToDir as FinalizeMergedRuns {
        input:
            files = [ annotated_bam, annotated_bai ],
            outdir = outdir + "/alignments"
    }
}
