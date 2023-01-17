version 1.0

import "tasks/bacteria/Bakta.wdl" as Bakta

workflow BaktaAnnotate {
    input {
        File bakta_db_tar
        File genome_fasta
        String? fname_prefix
    }

    call Bakta.BaktaAnnotate as Annotate {
        input:
            bakta_db_tar=bakta_db_tar,
            genome_fasta=genome_fasta,
            fname_prefix=fname_prefix
    }

    output {
        File tsv = Annotate.tsv
        File gff = Annotate.gff
        File genbank = Annotate.genbank
        File embl = Annotate.embl
        File ffn = Annotate.ffn
        File faa = Annotate.faa
        File hypotheticals_tsv = Annotate.hypotheticals_tsv
        File hypotheticals_faa = Annotate.hypotheticals_faa

        File summary = Annotate.summary
        File plot_png = Annotate.plot_png
        File plot_svg = Annotate.plot_svg
    }
}
