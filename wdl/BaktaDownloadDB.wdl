version 1.0

import "tasks/bacteria/Bakta.wdl" as Bakta
import "tasks/Finalize.wdl" as Finalize

workflow BaktaDownloadDB {
    input {
        String gcs_output_fname
    }

    String name = basename(gcs_output_fname)
    String outdir = basename(gcs_output_fname, name)

    call Bakta.BaktaDBDownload as Download { }
    call Finalize.FinalizeToFile as GCS {
        input:
            file = Download.bakta_db,
            outdir = outdir,
    }
}
