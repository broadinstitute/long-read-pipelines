version 1.0

##########################################################################################
# Top level workflow runner for Guppy.wdl, see there for more documentation 
##########################################################################################

import "tasks/Guppy.wdl" as Guppy
import "tasks/Finalize.wdl" as FF

workflow GuppyRunner {
    input {
        String gcs_fast5_dir
        String config = "dna_r9.4.1_450bps_hac.cfg"
        Boolean merge_fastqs = false
        String gcs_output_dir
    }

    call Guppy.Guppy {
        input:
            gcs_fast5_dir = gcs_fast5_dir,
            config        = config,
            merge_fastqs  = merge_fastqs
    }

    call FF.FinalizeToDir {
        input:
            files = Guppy.output_files,
            outdir = gcs_output_dir
    }

    output {
        Array[File] fastqs = Guppy.output_files
    }
}

