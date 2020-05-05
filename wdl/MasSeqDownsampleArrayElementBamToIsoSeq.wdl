version 1.0

import "tasks/TranscriptAnalysis/Preprocessing_Tasks.wdl" as TX_PRE
import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as UTILS
import "tasks/Finalize.wdl" as FF

workflow MasSeqDownsampleArrayElementBamToIsoSeq {

    meta {
        description :  "Downsample a given MAS-seq array element bam file into one containing only 1 read per ZMW (equivalent to IsoSeq)."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File array_element_bam
        String sample_name

        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/MasSeqDownsampleArrayElementBamToIsoSeq"
    }

    parameter_meta {
        array_element_bam : "Bam file containing MAS-seq array elements."
        sample_name : "Name of the sample being processed."
        gcs_out_root_dir  : "GCS bucket to store the output."
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    # Downsample by various factors:
    call TX_PRE.DownsampleToIsoSeqEquivalent as t01_DownsampleToIsoSeqEquivalent {
        input:
            array_element_bam = array_element_bam,
            prefix = sample_name + "_downsampled_to_isoseq"
    }

    ##########
    # store the results into designated bucket
    ##########

    String base_out_dir = outdir + "/" + sample_name + "/"

    call FF.FinalizeToDir as t09_FinalizeDownsampledToIsoseqReads {
        input:
            files = [t01_DownsampleToIsoSeqEquivalent.downsampled_bam],
            outdir = base_out_dir + "downsampled_to_isoseq"
    }
}
