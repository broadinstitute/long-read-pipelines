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
        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs."
    }

    input {
        File joint_vcf
        File joint_vcf_tbi
        String prefix

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ConvertToZarrStore/~{prefix}"

    # Gather across multiple input gVCFs
    call SGKit.ConvertToZarrStore as PerformZarrStoreConversion {
        input:
            vcf = joint_vcf,
            tbi = joint_vcf_tbi,
            prefix = prefix,
            outdir = outdir
    }

    ##########
    # store the results into designated bucket
    ##########

    output {
        String joint_zarr = PerformZarrStoreConversion.gcs_path
    }
}
