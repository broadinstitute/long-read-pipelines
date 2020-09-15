version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.34/wdl/tasks/Pacasus.wdl" as Pacasus
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.34/wdl/tasks/Finalize.wdl" as FF

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

