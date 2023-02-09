version 1.0

import "tasks/bacteria/Panaroo.wdl" as Panaroo

workflow RunPanaroo {
    input {
        Array[File] input_gffs
    }

    call Panaroo.Panaroo {
        input:
            input_gffs = input_gffs
    }

    output {
        File summary_stats = Panaroo.summary_stats
        File final_graph = Panaroo.final_graph
        File pre_filt_graph = Panaroo.pre_filt_graph
        File pangenome_reference = Panaroo.pangenome_reference

        File gene_data = Panaroo.gene_data
        File gene_presence_absence_rtab = Panaroo.gene_presence_absence_rtab
        File gene_presence_absence_csv = Panaroo.gene_presence_absence_csv
        File gene_presence_absence_roary = Panaroo.gene_presence_absence_roary

        File combined_DNA_CDS = Panaroo.combined_DNA_CDS
        File combined_protein_CDS = Panaroo.combined_protein_CDS
        File combined_protein_cdhit_out = Panaroo.combined_protein_cdhit_out
        File combined_protein_cdhit_out_cluster = Panaroo.combined_protein_cdhit_out_cluster

        File sv_triplets_presence_absence = Panaroo.sv_triplets_presence_absence
    }
}
