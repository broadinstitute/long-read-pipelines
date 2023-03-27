version 1.0

############################################################################################
## A workflow that performs joint calling on gVCFs (usually from DeepVariant) using GLNexus.
############################################################################################

import "../../../tasks/VariantCalling/GLNexus.wdl" as GLNexus
import "../../../tasks/Utility/Hail.wdl" as Hail
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow RunConvertToHailMT {
    input {
        File joint_gvcf
        File joint_gvcf_tbi
        String prefix

        String gcs_out_root_dir
    }

    parameter_meta {
        joint_gvcf:       "joint-called gVCF file"
        joint_gvcf_tbi:   ".tbi index for joint-called gVCF file"
        prefix:           "prefix for output Hail MatrixTable"
        gcs_out_root_dir: "GCS bucket in which to store the Hail MatrixTable"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/JointCallGVCFs/~{prefix}"

    # Gather across multiple input gVCFs
    call Hail.ConvertToHailMT {
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
        String joint_mt = ConvertToHailMT.gcs_path
    }
}
