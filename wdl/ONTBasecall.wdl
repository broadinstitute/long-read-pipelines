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

    call FF.FinalizeToDir as FinalizeFastqs {
        input:
            files = Guppy.output_files,
            outdir = gcs_out_root_dir
    }

    call FF.FinalizeToFile as FinalizeSequencingSummary {
        input:
            file = Guppy.sequencing_summary,
            outfile = gcs_out_root_dir + "/sequencing_summary.txt"
    }

    call FF.FinalizeToFile as FinalizeFinalSummary {
        input:
            file = Guppy.final_summary,
            outfile = gcs_out_root_dir + "/final_summary.txt"
    }

    output {
        String gcs_basecall_dir = FinalizeFastqs.gcs_dir
        File sequencing_summary = FinalizeSequencingSummary.gcs_path
        File final_summary = FinalizeFinalSummary.gcs_path
    }
}

