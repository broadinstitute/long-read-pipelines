version 1.0

import "tasks/bacteria/Bakta.wdl" as Bakta

workflow BaktaAnnotateBatch {
    input {
        File bakta_db_tar
        Int num_workers = 64
        String gcs_output_dir

        Array[File] genome_fastas
        Array[String] fname_prefixes
    }

    Int total = length(genome_fastas)
    Int batch_size = ceil(total / num_workers)

    String output_dir = sub(gcs_output_dir, "/+$", "")

    scatter(i in range(num_workers)) {
        call Bakta.BaktaAnnotateBatch as Annotate {
            input:
                bakta_db_tar=bakta_db_tar,
                output_dir=output_dir,
                all_genome_fastas=genome_fastas,
                all_fname_prefixes=fname_prefixes,

                worker=i,
                batch_size=batch_size
        }

    }

    output {
        Array[String] tsv = flatten(Annotate.tsv)
        Array[String] json = flatten(Annotate.json)
        Array[String] gff = flatten(Annotate.gff)
        Array[String] genbank = flatten(Annotate.genbank)
        Array[String] embl = flatten(Annotate.embl)
        Array[String] ffn = flatten(Annotate.ffn)
        Array[String] faa = flatten(Annotate.faa)
        Array[String] hypotheticals_tsv = flatten(Annotate.hypotheticals_tsv)
        Array[String] hypotheticals_faa = flatten(Annotate.hypotheticals_faa)

        Array[String] summary = flatten(Annotate.summary)
        Array[String] log = flatten(Annotate.log)
        Array[String] plot_png = flatten(Annotate.plot_png)
        Array[String] plot_svg = flatten(Annotate.plot_svg)
    }
}
