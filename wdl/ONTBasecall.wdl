version 1.0

import "tasks/Guppy.wdl" as Guppy
import "tasks/Finalize.wdl" as FF

workflow ONTBasecall {
    input {
        String gcs_fast5_dir
        String config = "dna_r9.4.1_450bps_hac.cfg"
        Array[String] barcode_kits = []
        String gcs_out_root_dir
    }

    call Guppy.Guppy {
        input:
            gcs_fast5_dir    = gcs_fast5_dir,
            config           = config,
            barcode_kits     = barcode_kits,
            gcs_out_root_dir = gcs_out_root_dir
    }

    output {
        String gcs_dir = Guppy.gcs_dir
        Array[File] sequencing_summaries = Guppy.sequencing_summaries
        Array[File] final_summaries = Guppy.final_summaries
        Array[String] barcodes = Guppy.barcodes
    }
}

