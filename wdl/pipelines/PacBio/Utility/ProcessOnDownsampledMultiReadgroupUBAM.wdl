version 1.0

import "ProcessOnInstrumentDemuxedChunk.wdl" as TreatPerReadGroup

import "../../../tasks/Utility/BAMutils.wdl" as BU
import "../../../tasks/Utility/GeneralUtils.wdl" as GU

workflow ProcessOnDownsampledMultiReadgroupUBAM {

    meta {
        description: "Align a downsampled multi-readgroup uBAM, perform QC check, and collect a few metrics."
    }

    parameter_meta {

        #########
        # inputs
        bam_SM_field:
        "value to place in the SM field of the resulting BAM header's @RG line"

        gcs_out_root_dir:
        "output files will be copied over there"

        qc_metrics_config_json:
        "A config json to for running the QC and metrics-collection sub-workflow 'AlignedBamQCandMetrics'"

        fingerprint_sample_id:
        "For fingerprint verification: the ID of the sample supposedly this BAM belongs to; note that the fingerprint VCF is assumed to be located at {fingerprint_store}/{fingerprint_sample_id}*.vcf(.gz)?"

        expected_sex_type:
        "If provided, triggers sex concordance check. Accepted value: [M, F, NA, na]"

        #########
        # outputs
    }

    input {
        String gcs_out_root_dir

        File  multiRGuBAM

        String bam_SM_field

        # args for optional QC subworkflows
        File? qc_metrics_config_json
        String? fingerprint_sample_id
        String? expected_sex_type

        File ref_map_file
        String disk_type = "SSD"
    }

    output {
        String last_processing_date = today.yyyy_mm_dd

        Array[File] readgroup_aligned_bams = eachRG.aligned_bam
        Array[File] readgroup_aligned_bais = eachRG.aligned_bai
        Array[File] readgroup_aligned_pbis = eachRG.aligned_pbi
        Array[String] readgroup_meta_movies = eachRG.movie
        Array[Float]  readgroup_metrics_contam = select_all(eachRG.contamination_est)

        Array[Map[String, String]?] readgroup_metrics_optional_fingerprint_check          = eachRG.fingerprint_check
        Array[Map[String, String]?] readgroup_metrics_optional_inferred_sex_info          = eachRG.inferred_sex_info
        Array[Map[String, String]?] readgroup_metrics_optional_methyl_tag_simple_stats    = eachRG.methyl_tag_simple_stats
        Array[Map[String, String]]  readgroup_metrics_optional_aBAM_metrics_files         = eachRG.aBAM_metrics_files
    }

    # split by read group
    call BU.SplitByRG { input:
        bam = multiRGuBAM,
        out_prefix = bam_SM_field,
        retain_rgless_records = false,
        sort_and_index = false
    }

    # align by each read group
    Array[Pair[String, File]] readgroup_2_bam = zip(SplitByRG.rg_ids, SplitByRG.split_bam)
    scatter (pair in readgroup_2_bam) {
        call TreatPerReadGroup.ProcessOnInstrumentDemuxedChunk as eachRG { input:
            gcs_out_root_dir = gcs_out_root_dir,

            uBAM = pair.right,

            readgroup_id = pair.left,
            bam_SM_field = bam_SM_field,

            platform = 'Revio',  # non-critical lie, this parameter was used for setting bam size lower bound threshold for failing QC

            # args for optional QC subworkflows
            qc_metrics_config_json = qc_metrics_config_json,
            fingerprint_sample_id = fingerprint_sample_id,
            expected_sex_type = expected_sex_type,

            ref_map_file = ref_map_file,
            disk_type = disk_type
        }
    }

    call GU.GetTodayDate as today {}
}
