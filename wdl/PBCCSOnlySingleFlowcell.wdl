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

workflow PBCCSOnlySingleFlowcell {
    input {
        String raw_reads_gcs_bucket

        String? sample_name
        Int num_shards = 300
        Boolean extract_uncorrected_reads = false

        String gcs_out_root_dir
    }

    parameter_meta {
        raw_reads_gcs_bucket:      "GCS bucket holding subreads BAMs (and other related files) holding the sequences to be CCS-ed"
        sample_name:               "[optional] name of sample this FC is sequencing"
        num_shards:                "[default-valued] number of sharded BAMs to create (tune for performance)"
        extract_uncorrected_reads: "[default-valued] extract reads that were not CCS-corrected to a separate file"
        gcs_out_root_dir :         "GCS bucket to store the corrected/uncorrected reads and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = raw_reads_gcs_bucket }

    # double scatter: one FC may generate multiple raw BAMs, we perform another layer scatter on each of these BAMs
    scatter (subread_bam in FindBams.subread_bams) {
        call PB.GetRunInfo { input: subread_bam = subread_bam }

        String SM  = select_first([sample_name, GetRunInfo.run_info["SM"]])
        String PU  = GetRunInfo.run_info["PU"]
        String ID  = PU
        String DIR = SM + "." + ID

        # break one raw BAM into fixed number of shards
        File subread_pbi = sub(subread_bam, ".bam$", ".bam.pbi")
        call PB.ShardLongReads { input: unaligned_bam = subread_bam, unaligned_pbi = subread_pbi, num_shards = num_shards }

        # then perform correction on each of the shard
        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads }

            if (extract_uncorrected_reads) {
                call PB.ExtractUncorrectedReads { input: subreads = subreads, consensus = CCS.consensus }
            }
        }

        # merge the corrected per-shard BAM/report into one, corresponding to one raw input BAM
        call Utils.MergeBams as MergeChunks { input: bams = CCS.consensus, prefix = "~{SM}.~{ID}" }
        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report }

        if (length(select_all(ExtractUncorrectedReads.uncorrected)) > 0) {
            call Utils.MergeBams as MergeUncorrectedChunks {
                input:
                    bams = select_all(ExtractUncorrectedReads.uncorrected),
                    prefix = "~{SM}.~{ID}.uncorrected"
            }
        }
    }

    # gather across (potential multiple) input raw BAMs
    if (length(FindBams.subread_bams) > 1) {
        call Utils.MergeBams as MergeRuns { input: bams = MergeChunks.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }
        call PB.MergeCCSReports as MergeAllCCSReports { input: reports = MergeCCSReports.report }

        if (length(select_all(MergeUncorrectedChunks.merged_bam)) > 0) {
            call Utils.MergeBams as MergeAllUncorrectedChunks {
                input:
                    bams = select_all(MergeUncorrectedChunks.merged_bam),
                    prefix = "~{SM}.~{ID}.uncorrected"
            }
        }
    }

    File ccs_bam = select_first([ MergeRuns.merged_bam, MergeChunks.merged_bam[0] ])
    File ccs_report = select_first([ MergeAllCCSReports.report, MergeCCSReports.report[0] ])

    if (extract_uncorrected_reads) {
        File? uncorrected_bam = select_first([ MergeAllUncorrectedChunks.merged_bam, MergeUncorrectedChunks.merged_bam[0] ])
    }

    ##########
    # store the results into designated bucket
    ##########

    call FF.FinalizeToDir as FinalizeCCSMetrics {
        input:
            files = [ ccs_report ],
            outdir = outdir + "/" + DIR[0] + "/metrics/ccs_metrics"
    }

    call FF.FinalizeToDir as FinalizeMergedRuns {
        input:
            files = [ ccs_bam ],
            outdir = outdir + "/" + DIR[0] + "/reads"
    }

    if (extract_uncorrected_reads) {
        call FF.FinalizeToDir as FinalizeMergedUncorrectedRuns {
            input:
                files = select_all([ uncorrected_bam ]),
                outdir = outdir + "/" + DIR[0] + "/reads"
        }
    }
}
