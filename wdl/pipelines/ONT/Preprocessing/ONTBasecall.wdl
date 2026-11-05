version 1.0

import "../../../tasks/Preprocessing/Guppy.wdl" as Guppy
import "../../../tasks/Utility/TerraUtils.wdl" as Terra

workflow ONTBasecall {

    meta {
        description: "Basecall ONT reads"
    }
    parameter_meta {
        gcs_fast5_dir: "GCS path to the directory containing fast5 files"
        config: "Guppy config file"
        barcode_kit: "Guppy barcode kit"
        gcs_out_root_dir: "Optional. GCS path to the root directory for output. When omitted, defaults to this run's Terra workspace bucket (auto-detected via GetWorkspaceMetadata), so no bucket need be supplied."
        prefix: "Prefix for output directory"
    }

    input {
        String gcs_fast5_dir
        String config = "dna_r10.4.1_e8.2_400bps_sup.cfg"
        String? barcode_kit
        String? gcs_out_root_dir
        String prefix
    }

    # Resolve the output root: caller-provided path, else auto-detect this run's Terra
    # workspace bucket (GetWorkspaceMetadata) so a bucket need not be supplied.
    if (!defined(gcs_out_root_dir)) {
        call Terra.GetWorkspaceMetadata as WorkspaceMetadata { input: }
    }
    String out_root = select_first([gcs_out_root_dir, WorkspaceMetadata.bucket])

    String outdir = sub(out_root, "/$", "") + "/ONTBasecall/~{prefix}"

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

