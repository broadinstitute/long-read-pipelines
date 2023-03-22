version 1.0

import "RunPanaroo.wdl" as Panaroo
import "IllumPlasmidSPAdes.wdl" as plasmidSPAdes
import "tasks/bacteria/Fusilli.wdl" as Fusilli

workflow IllumEvaluateFusilliPlasmids {
    input {
        File input_manifest
        File fusilli_db_tar
        File fusilli_db_index
        File panaroo_graph
        File panaroo_gene_data
        File panaroo_gene_presence_absence

        String sample_name
        File sample_fq1
        File sample_fq2
        Int graph_k = 49

        Array[String] true_plasmid_ids
        Array[String] true_plasmid_gff3s
    }

    call Fusilli.PrepareSample as PrepareSample {
        input:
            sample_name=sample_name,
            illumina_fq1=sample_fq1,
            illumina_fq2=sample_fq2,
            k=graph_k
    }

    call Fusilli.CallGenes as CallGenes {
        input:
            fusilli_db_tar=fusilli_db_tar,
            fusilli_db_index=fusilli_db_index,

            fusilli_prepared_sample=PrepareSample.fusilli_prepared_sample
    }

    Array[String] manifest_lines = read_lines(input_manifest)

    scatter(i in range(len(true_plasmid_ids))) {
        String tsv_row = "~{true_plasmid_ids[i]}\t~{true_plasmid_gff3s[i]}"
    }

    Array[String] manifest_with_truth = flatten([manifest_lines, tsv_row])

    # Create Panaroo graph with truth
    call Panaroo.RunPanaroo as PanarooTruth {
        input:
            input_manifest=write_lines(manifest_with_truth)
    }

    # Create Panaroo graph with plasmidSPAdes assembly
    call plasmidSPAdes.IllumPlasmidSPAdes as asmPlasmidSPAdes {
        input:
            illumina_fq1=sample_fq1,
            illumina_fq2=sample_fq2
    }

    Array[String] manifest_with_plasmidSPAdes = flatten([manifest_lines, ["plasmidSPAdes\t~{asmPlasmidSPAdes.contigs}"]])

    call Panaroo.RunPanaroo as PanarooSPAdes {
        input:
            input_manifest=write_lines(manifest_with_plasmidSPAdes)
    }

    output {
        File gene_predictions = CallGenes.gene_predictions

        File truth_panaroo_graph = PanarooTruth.final_graph
        File truth_panaroo_gene_data = PanarooTruth.gene_data
        File truth_panaroo_gene_presence_absence = PanarooTruth.gene_presence_absence_csv

        File spades_panaroo_graph = PanarooSPAdes.final_graph
        File spades_panaroo_gene_data = PanarooSPAdes.gene_data
        File spades_panaroo_gene_presence_absence = PanarooSPAdes.gene_presence_absence_csv
    }
}
