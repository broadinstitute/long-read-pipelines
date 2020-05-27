version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.12/wdl/tasks/Guppy.wdl" as Guppy
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.12/wdl/tasks/Finalize.wdl" as FF

workflow GuppyRunner {
    input {
        String gcs_fast5_dir
        String gcs_output_dir
    }

    call Guppy.Guppy {
        input:
            gcs_fast5_dir = gcs_fast5_dir
    }

    call FF.FinalizeToDir {
        input:
            files = Guppy.output_files,
            outdir = gcs_output_dir
    }

}

