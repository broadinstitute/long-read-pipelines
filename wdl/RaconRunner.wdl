version 1.0

##########################################################################################
# A task that polishes a draft assembly with long reads using Racon.
##########################################################################################

import "tasks/Racon.wdl" as Racon
import "tasks/Finalize.wdl" as FF

workflow RaconRunner {
    input {
        String reads
        String draft_assembly
        Int n_rounds

        String outdir
    }

    call Racon.RaconPolish {
        input:
            reads = reads,
            draft_assembly = draft_assembly,
            n_rounds = n_rounds
    }

    call FF.FinalizeToDir {
        input:
            files = RaconPolish.incremental_polished_assemblies,
            outdir = outdir
    }

}