version 1.0

import "../../../tasks/Preprocessing/GuppyGPU.wdl" as GuppyGPU
import "../../../tasks/Preprocessing/GuppyCPU.wdl" as GuppyCPU
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ONTBasecall {

    meta {
        description: "Basecall ONT reads"
    }
    parameter_meta {
        gcs_fast5_dir: "GCS path to the directory containing fast5 files"
        config: "Guppy config file"
        barcode_kit: "Guppy barcode kit"
        use_gpus: "If true, use the GPU version of the basecaller"
        gcs_out_root_dir: "GCS path to the root directory for output"
        prefix: "Prefix for output directory"
    }

    input {
        String gcs_fast5_dir
        String config = "dna_r10.4.1_e8.2_400bps_sup.cfg"
        String? barcode_kit
        Boolean use_gpus = true
        String gcs_out_root_dir
        String prefix
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTBasecall/~{prefix}"

    if (use_gpus) {
        call GuppyGPU.GuppyGPU {
            input:
                gcs_fast5_dir    = gcs_fast5_dir,
                config           = config,
                barcode_kit      = barcode_kit,
                gcs_out_root_dir = outdir
        }
    }

    if (!use_gpus) {
        call GuppyCPU.GuppyCPU {
            input:
                gcs_fast5_dir    = gcs_fast5_dir,
                config           = config,
                barcode_kit      = barcode_kit,
                gcs_out_root_dir = outdir
        }
    }

    output {
        String gcs_basecall_dir = select_first([GuppyGPU.gcs_dir, GuppyCPU.gcs_dir])
    }
}
