version 1.0

##########################################################################################
# A workflow that polishes a draft ONT assembly with Medaka.
##########################################################################################

import "tasks/Medaka.wdl" as Medaka
import "tasks/Finalize.wdl" as FF

workflow MedakaRunner {
    input {
        String basecalled_reads
        String draft_assembly
        String medaka_model

        String outdir
    }

    call Medaka.MedakaPolish {
        input:
            basecalled_reads = basecalled_reads,
            draft_assembly = draft_assembly,
            model = medaka_model
    }

    call FF.FinalizeToDir {
        input:
            files = [MedakaPolish.polished_assembly],
            outdir = outdir
    }

}