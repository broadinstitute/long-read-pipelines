## #

version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF

workflow PBCCS {

    meta {
        description: "A workflow that performs CCS correction on PacBio HiFi reads from a single flow cell. The workflow shards the subreads into clusters and performs CCS in parallel on each cluster. Ultimately, all the corrected reads (and uncorrected) are gathered into a single BAM. Various metrics are produced along the way."
    }

    input {
        Array[File] aligned_bams
        Array[File] aligned_bais

        String participant_name

        String gcs_out_root_dir
    }

    parameter_meta {
        aligned_bams:       "GCS path to aligned BAM files"
        aligned_bais:       "GCS path to aligned BAM file indices"

        participant_name:   "name of the participant from whom these samples were obtained"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCCS/~{participant_name}"

    # gather across (potential multiple) input CCS BAMs
    if (length(aligned_bams) > 1) {
        call Utils.MergeBams as MergeAllReads { input: bams = aligned_bams, prefix = participant_name }
    }

    File bam = select_first([MergeAllReads.merged_bam, aligned_bams[0]])
    File bai = select_first([MergeAllReads.merged_bai, aligned_bais[0]])

    call PB.PBIndex as IndexCCSUnalignedReads { input: bam = bam }
    File pbi = IndexCCSUnalignedReads.pbi

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
