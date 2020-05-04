version 1.0

import "tasks/Pacasus.wdl" as Pacasus
import "tasks/Finalize.wdl" as FF

workflow PacasusRunner {
    input {
        File input_fasta
        Int parallel_instances
        String output_dir
    }

    call Pacasus.Process {
        input:
            input_fasta = input_fasta,
            parallel_instances = parallel_instances
    }

    call FF.FinalizeToDir {
        input:
            files = [Process.processed_fasta],
            outdir = output_dir
    }
}