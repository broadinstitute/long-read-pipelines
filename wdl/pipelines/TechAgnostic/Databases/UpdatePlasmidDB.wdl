version 1.0

import "../../../tasks/Databases/PlasmidDB.wdl" as PlasmidDB

workflow UpdatePlasmidDB {
    input {
        String plsdb_fasta_url
        String plsdb_meta_url

        String gcs_output_dir
    }

    call PlasmidDB.UpdatePlasmidDB {
        input:
            plsdb_fasta_url=plsdb_fasta_url,
            plsdb_meta_url=plsdb_meta_url,
            gcs_output_directory=gcs_output_dir
    }
}
