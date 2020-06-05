version 1.0

import "tasks/Quast.wdl" as Quast
import "tasks/Finalize.wdl" as FF

workflow QuastRunner {
    input {
        File ref
        Array[File] assemblies

        String out_dir
    }

    call Quast.Quast {
        input:
            ref = ref,
            assemblies = assemblies
    }

    call FF.FinalizeToDir as FinalizeAssembly {
        input:
            files = [Quast.results],
            outdir = out_dir
    }
}
