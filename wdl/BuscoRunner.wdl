version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.36/wdl/tasks/Busco.wdl" as Busco
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.36/wdl/tasks/Finalize.wdl" as FF

workflow BuscoRunner {
    input {
        File assembly
        String lineage

        String out_dir
    }

    call Busco.Busco {
        input:
            assembly = assembly,
            lineage = lineage
    }

    call FF.FinalizeToDir as FinalizeAssembly {
        input:
            files = [Busco.output_tar],
            outdir = out_dir
    }
}
