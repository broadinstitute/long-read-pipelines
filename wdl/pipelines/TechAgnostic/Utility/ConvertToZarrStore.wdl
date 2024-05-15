version 1.0

import "../../../tasks/Utility/SGKit.wdl" as SGKit
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ConvertToZarrStore {

    meta {
        description: "A workflow that converts a joint-called VCF to a Zarr store."
    }

    parameter_meta {
        vcf: "Input joint called VCF file to convert to Zarr store format."
        joint_vcf_tbi: "VCF Index for `vcf`."
        prefix:           "prefix for output Zarr store"
        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs.  If no bucket is given, outputs will not be finalized and instead will remain in their native execution location."
    }

    input {
        File joint_vcf
        File joint_vcf_tbi
        String prefix

        String? gcs_out_root_dir
    }

    # Gather across multiple input gVCFs
    call SGKit.ConvertToZarrStore as PerformZarrStoreConversion {
        input:
            vcf = joint_vcf,
            tbi = joint_vcf_tbi,
            prefix = prefix
    }

    if (defined(gcs_out_root_dir)) {
        String concrete_gcs_out_root_dir = select_first([gcs_out_root_dir])
        call FF.FinalizeToFile as FinalizeOutputs {
            input:
                file = PerformZarrStoreConversion.zarr,
                outdir = concrete_gcs_out_root_dir
        }
    }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File joint_zarr = select_first([FinalizeOutputs.gcs_path, PerformZarrStoreConversion.zarr])
    }
}
