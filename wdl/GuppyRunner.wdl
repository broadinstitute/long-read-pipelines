version 1.0

##########################################################################################
# A workflow that runs the Guppy basecaller on ONT FAST5 files.
##########################################################################################

import "tasks/Guppy.wdl" as Guppy
import "tasks/Finalize.wdl" as FF

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

