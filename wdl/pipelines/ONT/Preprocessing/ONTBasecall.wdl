version 1.0

import "../../../tasks/Preprocessing/Guppy.wdl" as Guppy

workflow ONTBasecall {

    meta {
        description: "Basecall ONT reads"
    }
    parameter_meta {
        gcs_fast5_dir: "GCS path to the directory containing fast5 files"
        config: "Guppy config file"
        barcode_kit: "Guppy barcode kit"
        gcs_out_root_dir: "GCS path to the root directory for output"
        prefix: "Prefix for output directory"
    }

    input {
        String gcs_fast5_dir
        String config = "dna_r10.4.1_e8.2_400bps_sup.cfg"
        String? barcode_kit
        String gcs_out_root_dir
        String prefix
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTBasecall/~{prefix}"

    call Guppy.Guppy {
        input:
            gcs_fast5_dir    = gcs_fast5_dir,
            config           = config,
            barcode_kit      = barcode_kit,
            gcs_out_root_dir = outdir
    }

    output {
        String gcs_basecall_dir = Guppy.gcs_dir
    }
}

