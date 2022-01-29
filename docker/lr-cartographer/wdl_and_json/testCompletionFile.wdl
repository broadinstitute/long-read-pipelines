version 1.0

import "tasks/Finalize.wdl" as FF

workflow testCompletionFile {

    call FF.WriteCompletionFile {
        input:
            outdir = "gs://broad-dsde-methods-jonn/tmp/"
    }
}

