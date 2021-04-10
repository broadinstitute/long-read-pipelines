version 1.0

import "tasks/Guppy.wdl" as Guppy
import "tasks/Finalize.wdl" as FF

workflow ONTBasecall {
    input {
        String gcs_fast5_dir
        String config = "dna_r9.4.1_450bps_hac.cfg"
        String barcode_kit = "NA"
        String gcs_out_root_dir
        String prefix
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTBasecall/~{prefix}"

    call Guppy.Guppy {
        input:
            gcs_fast5_dir    = gcs_fast5_dir,
            config           = config,
            barcode_kit      = barcode_kit,
            gcs_out_root_dir = gcs_out_root_dir
    }

    output {
        String gcs_dir = Guppy.gcs_dir
        Array[File] sequencing_summaries = Guppy.sequencing_summaries
        Array[File] final_summaries = Guppy.final_summaries
        Array[String] barcodes = Guppy.barcodes
        Int num_fast5s = Guppy.num_fast5s
        Int num_pass_fastqs = Guppy.num_pass_fastqs
        Int num_fail_fastqs = Guppy.num_fail_fastqs
    }
}

