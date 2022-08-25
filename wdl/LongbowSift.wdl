version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF
import "tasks/Longbow.wdl" as LONGBOW

workflow LongbowSift {

    meta {
        description : "This workflow runs Longbow Sift on a given bam file."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File segmented_input_reads
        String sample_name

        String model = "mas_15_sc_10x5p_single_none"
        String validation_model = "10x_sc_10x5p_single_none"

        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/LongbowSift"
    }

    ##############################################################################################################

    call Utils.GetCurrentTimestampString as t_001_WdlExecutionStartTimestamp { input: }

    call LONGBOW.Sift as t_002_Sift {
        input:
            segmented_input_reads = segmented_input_reads,
            model = model,
            validation_model = validation_model,
            prefix = sample_name
    }

    ##############################################################################################################

    String outdir = sub(gcs_out_root_dir, "/$", "")
    String base_out_dir = outdir + "/" + sample_name + "/" + t_001_WdlExecutionStartTimestamp.timestamp_string

    # Finalize the final annotated, aligned array elements:
    call FF.FinalizeToDir as t_003_FinalizeQuantifiedArrayElements {
        input:
            files = [
                t_002_Sift.sifted_bam,
                t_002_Sift.sift_failed_bam,
                t_002_Sift.stats_tsv,
                t_002_Sift.summary_stats_tsv,
            ],
            outdir = base_out_dir
    }
}
