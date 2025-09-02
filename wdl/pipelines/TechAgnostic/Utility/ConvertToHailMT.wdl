version 1.0

import "../../../tasks/Utility/Hail.wdl" as Hail

workflow ConvertToHailMT {

    meta {
        description: "Convert a gVCF to a Hail MatrixTable"
    }
    parameter_meta {
        joint_gvcf:       "joint-called gVCF file"
        joint_gvcf_tbi:   ".tbi index for joint-called gVCF file"
        prefix:           "prefix for output Hail MatrixTable"
        gcs_out_root_dir: "GCS bucket in which to store the Hail MatrixTable"
    }

    input {
        File joint_gvcf
        File joint_gvcf_tbi
        String prefix

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/JointCallGVCFs/~{prefix}"

    # Gather across multiple input gVCFs
    call Hail.ConvertToHailMT as RunConvertToHailMT {
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
        String joint_mt = RunConvertToHailMT.mt_bucket_path
    }
}
