version 1.0

import "tasks/Structs.wdl"
import "tasks/bacteria/Bakta.wdl" as Bakta
import "tasks/Finalize.wdl" as Finalize

workflow BaktaDownloadDB {
    input {
        String gcs_output_fname
    }

    String name = basename(gcs_output_fname)
    String outdir = basename(gcs_output_fname, name)

    call Bakta.BaktaDBDownload as Download { }

    # Make sure finalize can handle the large database tar.bz2
    RuntimeAttr finalize_attr = {
        "disk_gb": 100
    }

    call Finalize.FinalizeToFile as GCS {
        input:
            file = Download.bakta_db,
            outdir = outdir,
            runtime_attr_override = finalize_attr
    }
}
