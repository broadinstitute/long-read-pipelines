version 1.0

import "../../../structs/Structs.wdl"
import "../../../tasks/GeneAnnotation/Bakta.wdl" as Bakta
import "../../../tasks/Utility/Finalize.wdl" as Finalize

workflow BaktaDownloadDB {
    input {
        String gcs_output_fname
    }

    String name = basename(gcs_output_fname)
    String outdir = sub(gcs_output_fname, "/([^/]+)$", "")

    call Bakta.BaktaDBDownload as Download {
        input:
            filename = name
    }

    call Finalize.FinalizeToFile as GCS {
        input:
            file = Download.bakta_db,
            outdir = outdir
    }
}
