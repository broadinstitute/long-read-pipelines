version 1.0

##########################################################################################
## A workflow that performs CCS correction on PacBio HiFi reads from a single flow cell.
## The workflow shards the subreads into clusters and performs CCS in parallel on each cluster.
## Ultimately, all the corrected reads (and uncorrected) are gathered into a single BAM.
## Various metrics are produced along the way.
##########################################################################################

import "../tasks/Utility/Utils.wdl" as Utils

workflow PBCCS {
    input {
        Array[File] aligned_bams

        String participant_name

        String gcs_out_root_dir
    }

    parameter_meta {
        aligned_bams:       "GCS path to aligned BAM files"

        participant_name:   "name of the participant from whom these samples were obtained"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCCS/~{participant_name}"

    # gather across (potential multiple) input CCS BAMs
    call Utils.MergeBams as MergeAllReads {
        input:
            bams = aligned_bams,
            outputBamName = "~{participant_name}.bam",
            outputBucket = outdir + "/alignments",
            pacBioBams = true
    }

    output {
        File aligned_bam = MergeAllReads.merged_bam
        File aligned_bai = MergeAllReads.merged_bai
        File aligned_pbi = select_first([MergeAllReads.merged_pbi])
    }
}
