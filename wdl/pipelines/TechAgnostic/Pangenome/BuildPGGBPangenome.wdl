version 1.0

import "../../../tasks/Pangenome/PGGB.wdl" as PGGB

workflow BuildPGGBPangenome {
    input {
        String pangenome_name

        Array[File] input_genomes
        String extra_pggb_args = ""
    }

    call PGGB.PrepareForPGGB as Prepare {
        input:
            output_fname_prefix=pangenome_name,
            input_genomes=input_genomes
    }

    call PGGB.PGGB as Build {
        input:
            pangenome_name=pangenome_name,
            pangenome_fasta_gz=Prepare.pangenome_fasta_gz,
            pangenome_fai=Prepare.pangenome_fai,
            extra_pggb_args=extra_pggb_args
    }

    output {
        File pangenome_fasta_gz = Prepare.pangenome_fasta_gz
        File pangenome_fai = Prepare.pangenome_fai

        File pangenome_gfa = Build.pangenome_gfa
        File pangenome_wfmash_paf = Build.pangenome_wfmash_paf

        File pangenome_odgi = Build.pangenome_odgi
        File pangenome_viz_1d = Build.pangenome_viz_1d

        File pangenome_odgi_layout = Build.pangenome_odgi_layout
        File pangenome_odgi_layout_tsv = Build.pangenome_odgi_layout_tsv
        File pangenome_viz_2d = Build.pangenome_viz_2d
    }
}