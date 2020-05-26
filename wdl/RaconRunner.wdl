version 1.0

import "tasks/Racon.wdl" as Racon
import "tasks/Finalize.wdl" as FF

workflow RaconRunner {
    input {
        String reads
        String draft_assembly

        String outdir
    }

    call Racon.RaconPolish {
        input:
            reads = reads,
            draft_assembly = draft_assembly
    }

    call FF.FinalizeToDir {
        input:
            files = [RaconPolish.polished_assembly],
            outdir = outdir
    }
    
}