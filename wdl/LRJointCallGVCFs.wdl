version 1.0

############################################################################################
## A workflow that performs joint calling on gVCFs (usually from DeepVariant) using GLNexus.
############################################################################################

import "tasks/GLNexus.wdl" as GLNexus
import "tasks/Finalize.wdl" as FF

workflow LRJointCallGVCFs {
    input {
        Array[File] gvcfs
        File? bed
        String prefix

        String gcs_out_root_dir
    }

    parameter_meta {
        gvcfs:            "GCS path to aligned BAM files"
        bed:              "three-column BED file with ranges to analyze"
        prefix:           "prefix for output joint-called gVCF and tabix index"
        gcs_out_root_dir: "GCS bucket to store the reads, variants, and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/JointCallGVCFs/~{prefix}"

    # Gather across multiple input gVCFs
    call GLNexus.JointCall { input: gvcfs = gvcfs, bed = bed, prefix = prefix }

    # Finalize
    call FF.FinalizeToFile as FinalizeGVCF { input: outdir = outdir, file = JointCall.joint_gvcf }
    call FF.FinalizeToFile as FinalizeTBI { input: outdir = outdir, file = JointCall.joint_gvcf_tbi }
    call FF.FinalizeToFile as FinalizeMT { input: outdir = outdir, file = JointCall.joint_mt_tar_gz }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File joint_gvcf = FinalizeGVCF.gcs_path
        File joint_gvcf_tbi = FinalizeTBI.gcs_path
        File joint_mt_tar_gz = FinalizeMT.gcs_path
    }
}
