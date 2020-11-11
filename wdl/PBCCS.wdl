version 1.0

##########################################################################################
## A workflow that performs CCS correction on PacBio HiFi reads from a single flow cell.
## The workflow shards the subreads into clusters and performs CCS in parallel on each cluster.
## Ultimately, all the corrected reads (and uncorrected) are gathered into a single BAM.
## Various metrics are produced along the way.
##########################################################################################

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF

workflow PBCCS {
    input {
        Array[File] bams

        String participant_name
        Int num_shards = 300
        Boolean extract_uncorrected_reads = false

        String? gcs_out_root_dir
    }

    parameter_meta {
        bams:                      "GCS path to raw subreads or CCS data"

        participant_name:          "name of the participant from whom these samples were obtained"
        num_shards:                "[default-valued] number of sharded BAMs to create (tune for performance)"
        extract_uncorrected_reads: "[default-valued] extract reads that were not CCS-corrected to a separate file"

        gcs_out_root_dir:          "[optional] GCS bucket to store the corrected/uncorrected reads, variants, and metrics files"
    }

    call Utils.GetDefaultDir { input: workflow_name = "PBCCS" }
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

            if (extract_uncorrected_reads) {
                call PB.ExtractUncorrectedReads { input: subreads = subreads, consensus = CCS.consensus }
            }
        }

        # merge the corrected per-shard BAM/report into one, corresponding to one raw input BAM
        call Utils.MergeBams as MergeCorrected { input: bams = CCS.consensus, prefix = "~{participant_name}.~{ID}.corrected" }

        if (length(select_all(ExtractUncorrectedReads.uncorrected)) > 0) {
            call Utils.MergeBams as MergeUncorrected {
                input:
                    bams = select_all(ExtractUncorrectedReads.uncorrected),
                    prefix = "~{participant_name}.~{ID}.uncorrected"
            }
        }

        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report, prefix = "~{participant_name}.~{ID}" }

        call FF.FinalizeToDir as FinalizeCCSReport {
            input:
                files = [ MergeCCSReports.report ],
                outdir = outdir + "/metrics/per_flowcell/" + ID + "/ccs_metrics"
        }
    }

    # gather across (potential multiple) input raw BAMs
    if (length(bams) > 1) {
        call Utils.MergeBams as MergeAllCorrected { input: bams = MergeCorrected.merged_bam, prefix = "~{participant_name}.corrected" }
        call PB.MergeCCSReports as MergeAllCCSReports { input: reports = MergeCCSReports.report }

        if (length(select_all(MergeUncorrected.merged_bam)) > 0) {
            call Utils.MergeBams as MergeAllUncorrected {
                input:
                    bams = select_all(MergeUncorrected.merged_bam),
                    prefix = "~{participant_name}.uncorrected"
            }
        }
    }

    File ccs_bam = select_first([ MergeAllCorrected.merged_bam, MergeCorrected.merged_bam[0] ])
    File ccs_report = select_first([ MergeAllCCSReports.report, MergeCCSReports.report[0] ])

    if (extract_uncorrected_reads) {
        File? uncorrected_bam = select_first([ MergeAllUncorrected.merged_bam, MergeUncorrected.merged_bam[0] ])
    }

    ##########
    # store the results into designated bucket
    ##########

    call FF.FinalizeToDir as FinalizeCCSMetrics {
        input:
            files = [ ccs_report ],
            outdir = outdir + "/metrics/combined/" + participant_name + "/ccs_metrics"
    }

    call FF.FinalizeToDir as FinalizeMergedRuns {
        input:
            files = [ ccs_bam ],
            outdir = outdir + "/reads"
    }

    if (extract_uncorrected_reads) {
        call FF.FinalizeToDir as FinalizeMergedUncorrectedRuns {
            input:
                files = select_all([ uncorrected_bam ]),
                outdir = outdir + "/reads"
        }
    }
}
