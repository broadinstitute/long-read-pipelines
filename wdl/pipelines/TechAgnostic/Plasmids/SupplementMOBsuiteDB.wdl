version 1.0

import "../../../tasks/Plasmids/MOBsuite.wdl" as MOBsuite

workflow SupplementMOBsuiteDB {
    input {
        File additional_plasmids_fasta
        File plasmid_host_tsv
    }

    call MOBsuite.CreateMOBsuiteDB as DB {
        input:
            additional_plasmids_fasta=additional_plasmids_fasta,
            plasmid_host_tsv=plasmid_host_tsv
    }

    output {
        File mobsuite_db = DB.db
    }
}