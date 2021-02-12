version 1.0

##########################################################################################
## A workflow that performs CCS correction on PacBio HiFi reads from a single sample,
## potentially involving multiple flowcells.  This workflow performs this activity without
## any data sharding whatsoever, but rather attempts to correct a single BAM file on a
## very large machine.  This is meant to mimic PacBio's default data processing on AWS, and
## is purposefully not optimized for GCP to serve as a point of comparison.
##
## DO NOT USE THIS WORKFLOW FOR PRODUCTION DATA PROCESSING.
##########################################################################################

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF

workflow PBCCSNoSharding {
    input {
        Array[File] bams

        String participant_name
        Int num_shards = 300
        Boolean extract_uncorrected_reads = false

        String gcs_out_root_dir
    }

    parameter_meta {
        bams:                      "GCS path to raw subreads or CCS data"

        participant_name:          "name of the participant from whom these samples were obtained"
        num_shards:                "[default-valued] number of sharded BAMs to create (tune for performance)"
        extract_uncorrected_reads: "[default-valued] extract reads that were not CCS-corrected to a separate file"

        gcs_out_root_dir:          "GCS bucket to store the corrected/uncorrected reads, variants, and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCCSNoSharding/" + participant_name

    # scatter over all sample BAMs
    scatter (bam in bams) {
        call PB.GetRunInfo { input: bam = bam }
        String ID = GetRunInfo.run_info["PU"]

        # perform correction on each BAM
        call PB.CCS {
            input:
                subreads = bam,
                cpus = 48,
                runtime_attr_override = {
                    'memory': 192,
                    'preemptible_tries': 0,
                    'max_retries': 0,
                }
        }

        call FF.FinalizeToDir as FinalizeCCSReport {
            input:
                files = [ CCS.report ],
                outdir = outdir + "/metrics/per_flowcell/" + ID + "/ccs_metrics"
        }
    }

    # gather across (potential multiple) input raw BAMs
    if (length(bams) > 1) {
        call Utils.MergeBams as MergeAllCorrected { input: bams = CCS.consensus, prefix = "~{participant_name}.corrected" }
        call PB.MergeCCSReports as MergeAllCCSReports { input: reports = CCS.report }
    }

    File ccs_bam = select_first([ MergeAllCorrected.merged_bam, CCS.consensus[0] ])
    File ccs_report = select_first([ MergeAllCCSReports.report, CCS.report[0] ])

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
}
