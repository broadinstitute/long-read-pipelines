version 1.0

import "tasks/Guppy.wdl" as Guppy
import "tasks/Finalize.wdl" as FF

workflow ONTBasecall {
    input {
        String gcs_fast5_dir
        String config = "dna_r9.4.1_450bps_hac.cfg"
        String gcs_out_root_dir
    }

    call Guppy.Guppy {
        input:
            gcs_fast5_dir = gcs_fast5_dir,
            config        = config
    }

    call FF.FinalizeToDir {
        input:
            files = Guppy.output_files,
            outdir = gcs_out_root_dir
    }

    output {
        String gcs_basecall_dir = FinalizeToDir.gcs_dir
    }
}

