version 1.0

import "../../../tasks/VariantCalling/GLNexus.wdl" as GLNexus
import "../../../tasks/Utility/Hail.wdl" as Hail
import "../../../tasks/Utility/VariantUtils.wdl" as VarUtils
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow LRConvertBCF {

    meta {
        description: "Convert a BCF file into a .vcf.gz file. Meant to temporarily handle some transient issues stemming from the LRJointCallGVCFs workflow. Should be removed eventually."
    }
    parameter_meta {
        joint_bcf: "The BCF file to convert"
        prefix: "The prefix to use for the output files"
        gcs_out_root_dir: "The root directory in GCS to store the output files"
    }

    input {
        File joint_bcf
        String prefix

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/JointCallGVCFs/~{prefix}"

    # Convert the joint .bcf callset into a single joint .vcf.bgz callset
    call VarUtils.ConcatVariants as ConcatBCFs { input: variant_files = [ joint_bcf ], prefix = prefix }

    call Hail.ConvertToHailMT {
        input:
            gvcf = ConcatBCFs.combined_vcf,
            tbi = ConcatBCFs.combined_vcf_tbi,
            prefix = prefix
    }

    # Finalize
    call FF.FinalizeToFile as FinalizeGVCF { input: outdir = outdir, file = ConcatBCFs.combined_vcf }
    call FF.FinalizeToFile as FinalizeTBI { input: outdir = outdir, file = ConcatBCFs.combined_vcf_tbi }
    call FF.FinalizeToFile as FinalizeHailMT { input: outdir = outdir, file = ConvertToHailMT.mt_tar }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File joint_vcf = FinalizeGVCF.gcs_path
        File joint_vcf_tbi = FinalizeTBI.gcs_path
        String joint_mt = FinalizeHailMT.gcs_path
    }
}
