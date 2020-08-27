version 1.0

import "tasks/Canu.wdl" as Canu
import "tasks/Quast.wdl" as Quast
import "tasks/Finalize.wdl" as FF

workflow MalariaAssembly {
    # guppy -> canu -> medaka -> quast?

    input {
        String gcs_fast5_dir
        String medaka_model

        String out_dir
    }

    call Guppy.Guppy {
        input:
            gcs_fast5_dir = gcs_fast5_dir
    }

    call Canu.CorrectTrimAssemble {
        input:
            output_file_prefix = "malaria",
            genome_size = "22.9m",
            reads = Guppy.output_files, # need to compress
            correct_corrected_error_rate = 0.15,
            trim_corrected_error_rate = 0.15,
            assemble_corrected_error_rate = 0.15
    }

    # add # of rounds
    call Medaka.MedakaPolish {
        input:
            basecalled_reads = Guppy.output_files, # compress
            draft_assembly = CorrectTrimAssemble.canu_contigs_fasta,
            model = medaka_model
    }

    call FF.FinalizeToDir {
        input:
            files = [MedakaPolish.polished_assembly],
            outdir = outdir
    }
}
