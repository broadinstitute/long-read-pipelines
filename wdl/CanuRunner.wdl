version 1.0

import "tasks/Canu.wdl" as Canu
import "tasks/Quast.wdl" as Quast
import "tasks/Finalize.wdl" as FF

workflow CanuRunner {
    input {
        File reads
        String output_file_prefix
        String genome_size
        Float correct_correctedErrorRate
        Float trim_correctedErrorRate
        Float assemble_correctedErrorRate

        String out_dir
    }

    call Canu.CorrectTrimAssemble {
        input:
            output_file_prefix = output_file_prefix,
            genome_size = genome_size,
            reads = reads,
            correct_corrected_error_rate = correct_correctedErrorRate,
            trim_corrected_error_rate = trim_correctedErrorRate,
            assemble_corrected_error_rate = assemble_correctedErrorRate
    }

    call FF.FinalizeToDir as FinalizeAssembly {
        input:
            files = [CorrectTrimAssemble.canu_contigs_fasta],
            outdir = out_dir
    }
}
