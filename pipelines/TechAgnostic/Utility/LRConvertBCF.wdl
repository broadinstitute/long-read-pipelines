version 1.0

############################################################################################
## A workflow that converts a BCF file into a .vcf.gz file. Meant to temporarily handle some
## transient issues stemming from the LRJointCallGVCFs workflow. Should be removed eventually.
############################################################################################

import "../../tasks/VariantCalling/GLNexus.wdl" as GLNexus
import "../../tasks/Utility/Hail.wdl" as Hail
import "../../tasks/Utility/Finalize.wdl" as FF

workflow LRConvertBCF {
    input {
        File joint_bcf
        String prefix

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/JointCallGVCFs/~{prefix}"

    # Convert the joint .bcf callset into a single joint .vcf.bgz callset
    call GLNexus.ConcatBCFs { input: bcfs = [ joint_bcf ], prefix = prefix }

    call Hail.ConvertToHailMT {
        input:
            gvcf = ConcatBCFs.joint_gvcf,
            tbi = ConcatBCFs.joint_gvcf_tbi,
            prefix = prefix,
            outdir = outdir
    }

    # Finalize
    call FF.FinalizeToFile as FinalizeGVCF { input: outdir = outdir, file = ConcatBCFs.joint_gvcf }
    call FF.FinalizeToFile as FinalizeTBI { input: outdir = outdir, file = ConcatBCFs.joint_gvcf_tbi }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File joint_gvcf = FinalizeGVCF.gcs_path
        File joint_gvcf_tbi = FinalizeTBI.gcs_path
        String joint_mt = ConvertToHailMT.gcs_path
    }
}
