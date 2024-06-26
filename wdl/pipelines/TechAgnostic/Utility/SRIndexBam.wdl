version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow SRIndexBam {

    meta {
        author: "Jonn Smith"
        description: "Index a given Bam file."
    }

    parameter_meta {
        bam: "Bam file for which to create an index."

        dir_prefix: "Directory prefix to use for finalized location."
        gcs_out_root_dir: "GCS Bucket into which to finalize outputs."
    }

    input {
        File bam

        String dir_prefix
        String gcs_out_root_dir
    }
    
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{dir_prefix}"

    call Utils.Index { input: bam = bam }
    call FF.FinalizeToFile as FinalizeBamIndex { input: outdir = outdir, file = Index.bai }

    output {
        File bai = FinalizeBamIndex.gcs_path
    }
}
