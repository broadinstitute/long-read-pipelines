version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.12/wdl/tasks/Quast.wdl" as Quast
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.12/wdl/tasks/Finalize.wdl" as FF

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
