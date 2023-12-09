version 1.0

import "../../../tasks/Utility/PBUtils.wdl" as PB
import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow PBFindAndSummarizeBAMs {

    meta {
        description: "A workflow that finds BAM files in a directory and emits some basic stats."
    }
    parameter_meta {
        gcs_input_dir:    "GCS path to input BAM files"
        length_threshold: "The minimum length of reads to consider in summary"

        gcs_out_root_dir: "GCS bucket to store the corrected/uncorrected reads, variants, and metrics files"
    }

    input {
        String gcs_input_dir
        Int length_threshold = 10000

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBFindAndSummarizeBAMs"

    call Utils.ListFilesOfType {
        input:
            gcs_dir = gcs_input_dir,
            suffixes = [".bam"],
            recurse = true
    }

    scatter (bam in ListFilesOfType.files[0]) {
        call PB.PBIndex { input: bam = bam }
        call PB.SummarizePBI { input: pbi = PBIndex.pbi, length_threshold = length_threshold }

        call FF.FinalizeToFile as FinalizePBISummary {
            input:
                outdir = outdir,
                file = SummarizePBI.results_file
        }
    }

    output {
    }
}