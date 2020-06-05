version 1.0

import "tasks/Pacasus.wdl" as Pacasus
import "tasks/Finalize.wdl" as FF

workflow PacasusRunner {
    input {
        File reads
        Int chunk_size_mb
        String output_dir
    }

    call Pacasus.Process {
        input:
            reads = reads,
            chunk_size_mb = chunk_size_mb
    }

    call FF.FinalizeToDir {
        input:
            files = [Process.processed_fasta],
            outdir = output_dir
    }
}

