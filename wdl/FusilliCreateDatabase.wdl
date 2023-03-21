version 1.0

import "tasks/bacteria/Fusilli.wdl" as Fusilli
import "RunPanaroo.wdl" as Panaroo

workflow FusilliCreateDatabase {
    input {
        String db_name
        File input_manifest

        Int index_k = 23
        Int? index_g
    }

    call Panaroo.RunPanaroo as PanarooSubWorkflow {
        input:
            input_manifest = input_manifest
    }

    call Fusilli.CreateDatabaseFromPanaroo as CreateDB {
        input:
            db_name = db_name,
            panaroo_graph = PanarooSubWorkflow.final_graph,
            panaroo_gene_data = PanarooSubWorkflow.gene_data,
            index_k=index_k,
            index_g=index_g
    }

    output {
        File summary_stats = PanarooSubWorkflow.summary_stats
        File final_graph = PanarooSubWorkflow.final_graph
        File pre_filt_graph = PanarooSubWorkflow.pre_filt_graph
        File pangenome_reference = PanarooSubWorkflow.pangenome_reference

        File gene_data = PanarooSubWorkflow.gene_data
        File gene_presence_absence_rtab = PanarooSubWorkflow.gene_presence_absence_rtab
        File gene_presence_absence_csv = PanarooSubWorkflow.gene_presence_absence_csv
        File gene_presence_absence_roary = PanarooSubWorkflow.gene_presence_absence_roary

        File combined_DNA_CDS = PanarooSubWorkflow.combined_DNA_CDS
        File combined_protein_CDS = PanarooSubWorkflow.combined_protein_CDS
        File combined_protein_cdhit_out = PanarooSubWorkflow.combined_protein_cdhit_out
        File combined_protein_cdhit_out_cluster = PanarooSubWorkflow.combined_protein_cdhit_out_cluster

        File sv_triplets_presence_absence = PanarooSubWorkflow.sv_triplets_presence_absence

        File fusilli_db_tar = CreateDB.fusilli_db_tar
        File fusilli_db_index = CreateDB.fusilli_db_index
    }

}
