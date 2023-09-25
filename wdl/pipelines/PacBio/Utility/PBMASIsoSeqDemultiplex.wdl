version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/PBUtils.wdl" as PB
import "../../../tasks/Preprocessing/Longbow.wdl" as Longbow
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow PBMASIsoSeqDemultiplex {

    meta {
        description: "Demultiplex a BAM file based on the UMI tag"
    }
    parameter_meta {
        bam:              "GCS path to BAM file"
        participant_name: "name of the participant from whom this sample was obtained"
        
        tag:              "[default valued] BAM tag on which to demultiplex (default: CB)"
        
        gcs_out_root_dir: "GCS bucket to store the demultiplexed BAMs"
    }

    input {
        File bam
        String participant_name

        String tag = "CB"

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBMASIsoSeqDemultiplex/~{participant_name}"

    # TODO: We should do this with 2 passes so that after we've ID'd the models involved, we go back and can refine the demux.  This is mostly for CLR reads.
    call Longbow.Demultiplex { input: bam = bam, tag = tag }

    scatter (demuxed_bam in Demultiplex.demuxed_bams) {
        call Utils.Index { input: bam = demuxed_bam }
        call PB.PBIndex { input: bam = demuxed_bam }

        # Finalize data
        call FF.FinalizeToFile as FinalizeBam { input: outdir = outdir, file = demuxed_bam, name = "~{participant_name}.bam" }
        call FF.FinalizeToFile as FinalizeBai { input: outdir = outdir, file = Index.bai,   name = "~{participant_name}.bam.bai" }
        call FF.FinalizeToFile as FinalizePbi { input: outdir = outdir, file = PBIndex.pbi, name = "~{participant_name}.bam.pbi" }
    }

    output {
        Array[File] demuxed_bams = FinalizeBam.gcs_path
        Array[File] demuxed_bais = FinalizeBai.gcs_path
        Array[File] demuxed_pbis = FinalizePbi.gcs_path
    }
}