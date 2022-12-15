version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF

workflow SRIndexBam {
    input {
        File bam
        String outdir
    }

    call Utils.Index { input: bam = bam }

    call FF.FinalizeToFile as FinalizeBamIndex { input: outdir = outdir, file = Index.bai }

    output {
        File bai = FinalizeBamIndex.gcs_path
    }
}