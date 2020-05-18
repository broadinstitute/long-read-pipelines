version 1.0

import "tasks/Pacasus.wdl" as Pacasus
import "tasks/Finalize.wdl" as FF

workflow PacasusSubRunner {
    input {
        File input_fasta
    }

    call Pacasus.RemovePalindromes {
        input:
            read_fasta = input_fasta
    }
}