version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/GCloudUtils.wdl" as GCPUtils
import "../../../tasks/Utility/VariantUtils.wdl" as VU

import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/Utility/ReadLengths.wdl" as ReLU

import "SmallVariantsBasicMetrics.wdl" as SingleSample

workflow JointVCFSmallVariantsBasicMetrics {
    meta {
        description:
        "Basic metrics on small variants joint VCF"
        note:
        "Since we aren't sure which QUAL filter will be good for you to look at, we run a range of potential QUAL filter values."
    }

    parameter_meta {
        pVCF:
        "Joint-call VCF."

        titrated_qualFilter_smallvar_stats:
        "Because we run a range of potential QUAL filter threshold values, the result is one entry per value, with the threshold value as the 'key'."
    }

    input {
        File pVCF
        File pVCFtbi

        String cohort_name

        String? hint_where_do_delocalize

        String gcs_out_root_dir
    }

    output {
        Array[Pair[String, Array[Pair[Int, Map[String, String]]]]] per_sample_titrated_qualFilter_smallvar_stats = sample_name_and_stats
    }

    String files_outdir = sub(gcs_out_root_dir, "/$", "") + "/JointVCFSmallVariantsBasicMetrics/~{cohort_name}"

    # here we take delocalization into our own hands, while at the same time ensuring call-caching
    # first, we need to make sure a hint about where to deloaclize is passed to the task that generates a large file or a large amount of files
    # then, the task, receiving the hint, performs the delocalization explicitly inside the task while side-stepping the usual delocalization mechanism by cromwell
    # lastly, to make sure call-caching is preserved, we must coerce the output type to File or Array[File]
    # call GCPUtils.GenerateDummyFile { }
    # call GCPUtils.GetParentDir { input: path = GenerateDummyFile.dummy }
    # String manual_delocalization_path_prefix = GetParentDir.parent_path + "/call-SplitJointCallVCFBySample"
    String? local_hint = files_outdir + "/" + hint_where_do_delocalize
    call VU.SplitJointCallVCFBySample { input:
        joint_vcf = pVCF,
        joint_vcf_tbi = pVCFtbi,
        hint_where_do_delocalize = local_hint # manual_delocalization_path_prefix # if the task is aliased, use call-<alias> here
    }
    if (defined(hint_where_do_delocalize)){
        call GCPUtils.ListFilesInDir as ListVCFs { input:
            path = select_first([local_hint]) + "/" + SplitJointCallVCFBySample.local_output_dir,
            pattern = "*.vcf.gz$"
        }
        call GCPUtils.ListFilesInDir as ListTBIs { input:
            path = select_first([local_hint]) + "/" + SplitJointCallVCFBySample.local_output_dir,
            pattern = "*.vcf.gz.tbi$"
        }
    }
    Array[File] splited_by_sample_vcfs = select_first([ListVCFs.files, SplitJointCallVCFBySample.vcfs])
    Array[File] splited_by_sample_tbis = select_first([ListTBIs.files, SplitJointCallVCFBySample.tbis])

    scatter (pair in zip(splited_by_sample_vcfs, splited_by_sample_tbis)) {
        File sample_vcf = pair.left
        File sample_tbi = pair.right
        call SingleSample.SmallVariantsBasicMetrics {
            input:
                vcf = sample_vcf,
                tbi = sample_tbi,
                gcs_out_root_dir = files_outdir
        }
        call VU.GetVCFSampleName { input: fingerprint_vcf = sample_vcf }
        String sample_name = GetVCFSampleName.sample_name
        Array[Pair[Int, Map[String, String]]] sample_stats = SmallVariantsBasicMetrics.titrated_qualFilter_smallvar_stats
    }
    Array[Pair[String, Array[Pair[Int, Map[String, String]]]]] sample_name_and_stats = zip(sample_name, sample_stats)
}
