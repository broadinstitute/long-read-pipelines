version 1.0

#################################################################
## A workflow that converts a joint-called gVCF to a Zarr store.
#################################################################

import "../../../tasks/Utility/SGKit.wdl" as SGKit
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ConvertToZarrStore {
    input {
        File joint_gvcf
        File joint_gvcf_tbi
        String prefix

        String gcs_out_root_dir
    }

    parameter_meta {
        joint_gvcf:       "joint-called gVCF file"
        joint_gvcf_tbi:   ".tbi index for joint-called gVCF file"
        prefix:           "prefix for output Zarr store"
        gcs_out_root_dir: "GCS bucket in which to store the Zarr store"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/JointCallGVCFs/~{prefix}"

    # Gather across multiple input gVCFs
    call SGKit.ConvertToZarrStore as PerformZarrStoreConversion {
        input:
            gvcf = joint_gvcf,
            tbi = joint_gvcf_tbi,
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
