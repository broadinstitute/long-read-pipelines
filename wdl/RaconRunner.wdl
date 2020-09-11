version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.30/wdl/tasks/Racon.wdl" as Racon
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.30/wdl/tasks/Finalize.wdl" as FF

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