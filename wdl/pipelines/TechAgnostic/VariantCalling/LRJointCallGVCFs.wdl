version 1.0

import "../../../tasks/VariantCalling/GLNexus.wdl" as GLNexus
import "../../../tasks/Utility/Hail.wdl" as Hail
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow LRJointCallGVCFs {

    meta {
        description: "A workflow that performs joint calling on gVCFs (usually from DeepVariant) using GLNexus."
    }
    parameter_meta {
        gvcfs:            "GCS paths to gVCF files"
        tbis:             "GCS paths to gVCF tbi files"
        ref_map_file:     "table indicating reference sequence and auxillary file locations"
        prefix:           "prefix for output joint-called gVCF and tabix index"
        gcs_out_root_dir: "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        Array[File] gvcfs
        Array[File] tbis
        File ref_map_file

        String prefix

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/JointCallGVCFs/~{prefix}"

    Map[String, String] ref_map = read_map(ref_map_file)

    # Gather across multiple input gVCFs
    call GLNexus.JointCall {
        input:
            gvcfs = gvcfs,
            tbis = tbis,
            dict = ref_map['dict'],
            prefix = prefix
    }

    call Hail.ConvertToHailMT {
        input:
            gvcf = JointCall.joint_gvcf,
            tbi = JointCall.joint_gvcf_tbi,
            prefix = prefix
    }

    # Finalize
    call FF.FinalizeToFile as FinalizeGVCF { input: outdir = outdir, file = JointCall.joint_gvcf }
    call FF.FinalizeToFile as FinalizeTBI { input: outdir = outdir, file = JointCall.joint_gvcf_tbi }
    call FF.FinalizeToFile as FinalizeMT { input: outdir = outdir, file = ConvertToHailMT.mt_tar }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File joint_gvcf = FinalizeGVCF.gcs_path
        File joint_gvcf_tbi = FinalizeTBI.gcs_path
        String joint_mt = FinalizeMT.gcs_path
    }
}
