version 1.0

import "../../../tasks/Pangenome/Panaroo.wdl" as Panaroo

workflow RunPanaroo {
    input {
        File input_manifest
    }

    call Panaroo.Panaroo as CallPanaroo {
        input:
            input_manifest = input_manifest
    }

    output {
        File summary_stats = CallPanaroo.summary_stats
        File final_graph = CallPanaroo.final_graph
        File pre_filt_graph = CallPanaroo.pre_filt_graph
        File pangenome_reference = CallPanaroo.pangenome_reference

        File gene_data = CallPanaroo.gene_data
        File gene_presence_absence_rtab = CallPanaroo.gene_presence_absence_rtab
        File gene_presence_absence_csv = CallPanaroo.gene_presence_absence_csv
        File gene_presence_absence_roary = CallPanaroo.gene_presence_absence_roary

        File combined_DNA_CDS = CallPanaroo.combined_DNA_CDS
        File combined_protein_CDS = CallPanaroo.combined_protein_CDS
        File combined_protein_cdhit_out = CallPanaroo.combined_protein_cdhit_out
        File combined_protein_cdhit_out_cluster = CallPanaroo.combined_protein_cdhit_out_cluster

        File sv_triplets_presence_absence = CallPanaroo.sv_triplets_presence_absence
    }
}
