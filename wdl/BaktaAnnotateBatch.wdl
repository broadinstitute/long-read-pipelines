version 1.0

import "tasks/bacteria/Bakta.wdl" as Bakta

workflow BaktaAnnotateBatch {
    input {
        File bakta_db_tar
        Int num_workers = 64
        String gcs_output_dir

        Array[String] plasmid_ids
        Array[File] genome_fastas
    }

    Int total = length(genome_fastas)
    Int batch_size = ceil(total / num_workers)

    String output_dir = sub(gcs_output_dir, "/+$", "")

    scatter(i in range(num_workers)) {
        call Bakta.BaktaAnnotateBatch as Annotate {
            input:
                bakta_db_tar=bakta_db_tar,
                output_dir=output_dir,
                plasmid_ids=plasmid_ids,
                all_genome_fastas=genome_fastas,

                worker=i,
                batch_size=batch_size
        }

    }

    call Bakta.CreateTerraDataTSV {
        input:
            plasmid_ids=plasmid_ids,
            tsv = flatten(Annotate.tsv),
            json = flatten(Annotate.json),
            gff = flatten(Annotate.gff),
            genbank = flatten(Annotate.genbank),
            embl = flatten(Annotate.embl),
            ffn = flatten(Annotate.ffn),
            faa = flatten(Annotate.faa),
            hypotheticals_tsv = flatten(Annotate.hypotheticals_tsv),
            hypotheticals_faa = flatten(Annotate.hypotheticals_faa),
            summary = flatten(Annotate.summary),
            log = flatten(Annotate.log),
            plot_png = flatten(Annotate.plot_png),
            plot_svg = flatten(Annotate.plot_svg),

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
