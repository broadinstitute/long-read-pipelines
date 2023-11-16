version 1.0

##########################################################################################
## A workflow that performs CCS correction on PacBio HiFi reads from a single flow cell.
## The workflow shards the subreads into clusters and performs CCS in parallel on each cluster.
## Ultimately, all the corrected reads (and uncorrected) are gathered into a single BAM.
## Various metrics are produced along the way.
##########################################################################################

import "../tasks/Utility/Utils.wdl" as Utils
import "../tasks/Utility/Finalize.wdl" as FF

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
            prefix = participant_name,
            pacBioBams = true
    }
    File bam = MergeAllReads.merged_bam
    File bai = MergeAllReads.merged_bai
    File pbi = select_first([MergeAllReads.merged_pbi])

    # Finalize
    String dir = outdir + "/alignments"

    call FF.FinalizeToFile as FinalizeBam { input: outdir = dir, file = bam, name = "~{participant_name}.bam" }
    call FF.FinalizeToFile as FinalizeBai { input: outdir = dir, file = bai, name = "~{participant_name}.bam.bai" }
    call FF.FinalizeToFile as FinalizePbi { input: outdir = dir, file = pbi, name = "~{participant_name}.bam.pbi" }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File aligned_bam = FinalizeBam.gcs_path
        File aligned_bai = FinalizeBai.gcs_path
        File aligned_pbi = FinalizePbi.gcs_path
    }
}
