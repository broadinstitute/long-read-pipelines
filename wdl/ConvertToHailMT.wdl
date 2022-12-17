version 1.0

############################################################################################
## A workflow that performs joint calling on gVCFs (usually from DeepVariant) using GLNexus.
############################################################################################

import "tasks/GLNexus.wdl" as GLNexus
import "tasks/Hail.wdl" as Hail
import "tasks/Finalize.wdl" as FF

workflow ConvertToHailMT {
    input {
        File joint_gvcf
        File joint_gvcf_tbi
        String prefix

        String gcs_out_root_dir
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
