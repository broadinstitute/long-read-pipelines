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
        background_sample_gvcfs: "Array of GVCFs to use as background samples for joint calling."
        background_sample_gvcf_indices: "Array of GVCF index files for `background_sample_gvcfs`.  Order should correspond to that in `background_sample_gvcfs`."
        gcs_out_root_dir: "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        Array[File] gvcfs
        Array[File] tbis
        File ref_map_file

        String prefix

        Array[Array[File]]? background_sample_gvcfs
        Array[Array[File]]? background_sample_gvcf_indices

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/JointCallGVCFs/~{prefix}"

    Map[String, String] ref_map = read_map(ref_map_file)

    # Gather across multiple input gVCFs
    call GLNexus.JointCall {
        input:
            gvcfs = gvcfs,
            tbis = tbis,
            background_sample_gvcfs = background_sample_gvcfs,
            background_sample_gvcf_indices = background_sample_gvcf_indices,
            dict = ref_map['dict'],
            force_add_missing_dp = true,
            prefix = prefix
    }

    call Hail.ConvertToHailMT {
        input:
            gvcf = JointCall.joint_vcf,
            tbi = JointCall.joint_vcf_tbi,
            prefix = prefix
    }

    # Finalize
    call FF.FinalizeToFile as FinalizeVCF { input: outdir = outdir, file = JointCall.joint_vcf }
    call FF.FinalizeToFile as FinalizeTBI { input: outdir = outdir, file = JointCall.joint_vcf_tbi }
    call FF.FinalizeToFile as FinalizeMT { input: outdir = outdir, file = ConvertToHailMT.mt_tar }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File joint_vcf = FinalizeVCF.gcs_path
        File joint_vcf_tbi = FinalizeTBI.gcs_path
        String joint_mt = FinalizeMT.gcs_path
    }
}
