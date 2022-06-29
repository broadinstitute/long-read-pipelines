version 1.0

import "tasks/Racon.wdl" as Racon



workflow mito_racon {

    input {
        File reads
        File draft_assembly
        Int n_rounds
    }

    parameter_meta {
        reads:          "long reads to polish the draft assembly with"
        draft_assembly: "draft to be polished"
        n_rounds: "Number of times to run Racon"
    }

    call Racon.RaconPolish as RaconPolish {
        input:
            reads = reads,
            draft_assembly = draft_assembly,
            n_rounds = n_rounds
    }

    output {
        File final_polished_assembly = RaconPolish.final_polished_assembly
        Array[File] incremental_polished_assemblies =  RaconPolish.incremental_polished_assemblies
    }
}