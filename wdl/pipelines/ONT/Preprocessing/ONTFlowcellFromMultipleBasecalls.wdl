version 1.0

import "../../../tasks/Utility/GeneralUtils.wdl" as GU

import "../../TechAgnostic/Utility/MergeSampleBamsAndCollectMetrics.wdl" as MERGE

workflow ONTFlowcellFromMultipleBasecalls {
    meta {
        description:
        "For merging aligned BAMs from the multiple basecall directories of a single flowcell (assumed to be from a single sample)."
    }
    input {
        Array[File] aligned_bams
        Array[File] aligned_bais
        Boolean bams_suspected_to_contain_dup_record = true

        String flowcell

        File? qc_metrics_config_json
        String? fingerprint_sample_id
        String? expected_sex_type

        File ref_map_file

        String gcs_out_root_dir
    }

    parameter_meta {
        bams_suspected_to_contain_dup_record: "if the multiple basecalls are suspected to have duplicates amongst them"
        ref_map_file:       "table indicating reference sequence and auxillary file locations"
        bed_to_compute_coverage: "A BED listing regions where each will get a coverage summary"
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTFlowcell/~{flowcell}"

    # yes, we are using the sample-merging sub-workflow here off-label, but the two really share the same logic
    call MERGE.Work as MergeAndMetrics {
        input:
            gcs_out_dir = outdir,

            sample_name = flowcell,  # this is basically the only place where the two differ
            aligned_bams = aligned_bams,
            aligned_bais = aligned_bais,

            ref_map_file = ref_map_file,

            tech = 'ONT',
            bams_suspected_to_contain_dup_record = bams_suspected_to_contain_dup_record,

            qc_metrics_config_json = qc_metrics_config_json,
            fingerprint_sample_id = fingerprint_sample_id,
            expected_sex_type = expected_sex_type
    }
    # todo: aggregate ONT-specific raw reads stats from basecall directories

    call GU.GetTodayDate as today {}

    output {
        String last_process_date = today.yyyy_mm_dd

        ########################################
        File aligned_bam = MergeAndMetrics.aligned_bam
        File aligned_bai = MergeAndMetrics.aligned_bai

        Float coverage = MergeAndMetrics.coverage

        ########################################
        # QC/metrics
        Map[String, Float] nanoplot_summ             = MergeAndMetrics.nanoplot_summ
        Map[String, Float] sam_flag_stats            = MergeAndMetrics.sam_flag_stats

        # fingerprint
        Map[String, String]? fingerprint_check       = MergeAndMetrics.fingerprint_check

        # contam
        Float? contamination_est                     = MergeAndMetrics.contamination_est

        # sex concordance
        Map[String, String]? inferred_sex_info       = MergeAndMetrics.inferred_sex_info

        # methyl
        Map[String, String]? methyl_tag_simple_stats = MergeAndMetrics.methyl_tag_simple_stats

        # file-based QC/metrics outputs all packed into a finalization map
        Map[String, String] aBAM_metrics_files       = MergeAndMetrics.aBAM_metrics_files
    }
}