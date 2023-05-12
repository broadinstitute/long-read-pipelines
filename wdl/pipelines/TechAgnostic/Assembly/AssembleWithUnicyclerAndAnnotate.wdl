version 1.0

import "AssembleWithUnicycler.wdl" as Unicycler
import "../Annotation/BaktaAnnotateSingle.wdl" as Bakta

workflow AssmebleWithUnicyclerAndAnnotate {
    input {
        String sample_name

        File illumina_fq1
        File illumina_fq2
        File? long_reads

        File bakta_db_tar
    }

    call Unicycler.AssembleWithUnicycler as Assembly {
        input:
            sample_name=sample_name,
            illumina_fq1=illumina_fq1,
            illumina_fq2=illumina_fq2,
            long_reads=long_reads,
    }

    call Bakta.BaktaAnnotateSingle as Annotate {
        input:
            bakta_db_tar=bakta_db_tar,
            genome_fasta=Assembly.assembly_fasta,
            fname_prefix=sample_name
    }

    output {
        File processed_fq1 = Assembly.processed_fq1
        File processed_fq2 = Assembly.processed_fq2
        File unpaired_fq = Assembly.unpaired_fq
        File fastp_report = Assembly.fastp_report
        File fastp_json = Assembly.fastp_json

        File assembly_fasta = Assembly.assembly_fasta
        File assembly_gfa = Assembly.assembly_gfa
        File unicycler_log = Assembly.unicycler_log
        Array[File] unicycler_intermediate_graphs = Assembly.unicycler_intermediate_graphs

        File annot_tsv = Annotate.tsv
        File annot_json = Annotate.json
        File annot_gff3 = Annotate.gff
        File annot_genbank = Annotate.genbank
        File annot_embl = Annotate.embl
        File cds_ffn = Annotate.ffn
        File protein_faa = Annotate.faa
        File hypotheticals_tsv = Annotate.hypotheticals_tsv
        File hypotheticals_faa = Annotate.hypotheticals_faa

        File annot_summary = Annotate.summary
        File annot_log = Annotate.log
        File plot_png = Annotate.plot_png
        File plot_svg = Annotate.plot_svg
    }
}
