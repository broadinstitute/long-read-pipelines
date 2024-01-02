version 1.0

import "../../../tasks/Utility/GeneralUtils.wdl" as GU

import "../../TechAgnostic/Utility/MergeSampleBamsAndCollectMetrics.wdl" as MERGE
import "../../TechAgnostic/VariantCalling/CallVariantsReadBased.wdl" as VAR

workflow PBCCSWholeGenome {

    meta {
        description: "A workflow that performs single sample variant calling on PacBio HiFi reads from one or more SMRT cells. The workflow merges multiple SMRT cells into a single BAM prior to variant calling."
    }
    parameter_meta {
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"

        aligned_bams:       "GCS path to aligned BAM files"
        aligned_bais:       "GCS path to aligned BAM file indices"
        sample_name:        "sample name as encoded in the bams"

        platform: "PacBio platform used for generating the data; accepted value: [Sequel, Revio]"

        ref_map_file: "table indicating reference sequence and auxillary file locations"
        ref_scatter_interval_list_locator: "A file holding paths to interval_list files, used for custom sharding the of the input BAM; when not provided, will shard WG by contig (possibly slower)"
        ref_scatter_interval_list_ids: "A file that gives short IDs to the interval_list files; when not provided, will shard WG by contig (possibly slower)"

        qc_metrics_config_json:
        "A config json to for running the QC and metrics-collection sub-workflow 'AlignedBamQCandMetrics'"

        fingerprint_sample_id:
        "For fingerprint verification: the ID of the sample supposedly this BAM belongs to; note that the fingerprint VCF is assumed to be located at {fingerprint_store}/{fingerprint_sample_id}*.vcf(.gz)?"

        expected_sex_type:
        "If provided, triggers sex concordance check. Accepted value: [M, F, NA, na]"

        call_svs:               "whether to call SVs"
        pbsv_discover_per_chr:  "Run the discover stage of PBSV per chromosome"

        call_small_variants: "whether to call small variants"

        run_clair3:  "to turn on Clair3 analysis or not (non-trivial increase in cost and runtime)"

        use_margin_for_tagging: "if false, will use margin-phased small-variant VCF for haplotagging the BAM; applicable only when input data isn't ONT data with pore older than R10.4"

        gcp_zones: "which Google Cloud Zone to use (this has implications on how many GPUs are available and egress costs, so configure carefully)"

        # metrics outputs
        nanoplot_summ:
        "Summary on alignment metrics provided by Nanoplot (todo: study the value of this output)"

        sam_flag_stats:
        "SAM flag stats"
        fingerprint_check:
        "Summary on (human) fingerprint checking results"
        contamination_est:
        "cross-(human)individual contamination estimation by VerifyBAMID2"
        inferred_sex_info:
        "Inferred sex concordance information if expected sex type is provided"
        methyl_tag_simple_stats:
        "Simple stats on the reads with & without SAM methylation tags (MM/ML)."
        aBAM_metrics_files:
        "A map where keys are summary-names and values are paths to files generated from the various QC/metrics tasks"

        # variants outputs
        dv_g_vcf: "DeepVariant gVCF; available for CCS data and ONT data generated with pores >= R10.4."
        dv_g_tbi: "Index for DeepVariant ; available for CCS data and ONT data generated with pores >= R10.4."
        dv_margin_phased_vcf: "Phased DeepVariant VCF genrated with Margin; available for CCS data and ONT data generated with pores >= R10.4."
        dv_margin_phased_tbi: "Index for phased DeepVariant VCF genrated with Margin; available for CCS data and ONT data generated with pores >= R10.4."
        dv_vcf_margin_phasing_stats_tsv: "Phasing stats (TSV format) of phased DeepVariant VCF genrated with Margin; available for CCS data and ONT data generated with pores >= R10.4."
        dv_vcf_margin_phasing_stats_gtf: "Phasing stats (GTF format) of phased DeepVariant VCF genrated with Margin; available for CCS data and ONT data generated with pores >= R10.4."
        dv_whatshap_phased_vcf: "Phased DeepVariant VCF genrated with WhatsHap; available for CCS data and ONT data generated with pores >= R10.4."
        dv_whatshap_phased_tbi: "Index for phased DeepVariant VCF genrated with WhatsHap; available for CCS data and ONT data generated with pores >= R10.4."
        dv_vcf_whatshap_phasing_stats_tsv: "Phasing stats (TSV format) of phased DeepVariant VCF genrated with WhatsHap; available for CCS data and ONT data generated with pores >= R10.4."
        dv_vcf_whatshap_phasing_stats_gtf: "Phasing stats (GTF format) of phased DeepVariant VCF genrated with WhatsHap; available for CCS data and ONT data generated with pores >= R10.4."

        dv_nongpu_resources_usage_visual: "Resource usage monitoring log visualization for DV (per shard); available for CCS data and ONT data generated with pores >= R10.4."
    }

    input {
        String gcs_out_root_dir

        # sample specific
        String sample_name
        Array[File] aligned_bams
        Array[File] aligned_bais
        String platform

        # reference-specific
        File ref_map_file
        File? ref_scatter_interval_list_locator
        File? ref_scatter_interval_list_ids

        # variant-calling user choice
        Boolean call_svs = true
        Boolean pbsv_discover_per_chr = true
        Int minsvlen = 50

        Boolean call_small_variants = true
        Boolean run_clair3 = false
        Boolean use_margin_for_tagging = true
        Int dv_threads = 16
        Int dv_memory = 40
        Boolean use_gpu = false

        # for QC/metrics
        File? qc_metrics_config_json
        String? fingerprint_sample_id
        String? expected_sex_type

        Array[String] gcp_zones = ['us-central1-a', 'us-central1-b', 'us-central1-c', 'us-central1-f']
    }

    String workflow_name = "PBCCSWholeGenome"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_name}/~{sample_name}"

    ###########################################################
    call MERGE.Work as MergeAndMetrics {
        input:
            gcs_out_dir = outdir,

            sample_name = sample_name,
            aligned_bams = aligned_bams,
            aligned_bais = aligned_bais,

            ref_map_file = ref_map_file,

            tech = platform,
            bams_suspected_to_contain_dup_record = false,

            qc_metrics_config_json = qc_metrics_config_json,
            fingerprint_sample_id = fingerprint_sample_id,
            expected_sex_type = expected_sex_type,

            run_seqkit_stats = false
    }

    ###########################################################
    if (call_svs || call_small_variants) {
        call VAR.CallVariants {
            input:
                gcs_out_dir = outdir,
                bam    = MergeAndMetrics.aligned_bam,
                bai    = MergeAndMetrics.aligned_bai,
                prefix = sample_name,

                is_ont = false,
                is_r10_4_pore_or_later = false,
                model_for_dv_andor_pepper = 'PACBIO',

                ref_map_file = ref_map_file,
                ref_scatter_interval_list_locator = ref_scatter_interval_list_locator,
                ref_scatter_interval_list_ids = ref_scatter_interval_list_ids,

                call_svs = call_svs,
                pbsv_discover_per_chr = pbsv_discover_per_chr,
                minsvlen = minsvlen,

                call_small_variants = call_small_variants,
                run_clair3 = run_clair3,
                use_margin_for_tagging = use_margin_for_tagging,
                dv_threads = dv_threads,
                dv_memory = dv_memory,
                use_gpu = use_gpu,

                gcp_zones = gcp_zones
        }
    }

    ###########################################################
    call GU.GetTodayDate as today {}

    ###########################################################
    output {
        String last_processing_date = today.yyyy_mm_dd

        ########################################
        File aligned_bam = MergeAndMetrics.aligned_bam
        File aligned_bai = MergeAndMetrics.aligned_bai
        File aligned_pbi = select_first([MergeAndMetrics.aligned_pbi])

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
        Map[String, String] aBAM_metrics_files      = MergeAndMetrics.aBAM_metrics_files

        ########################################
        # variants
        File? pbsv_vcf = CallVariants.pbsv_vcf
        File? pbsv_tbi = CallVariants.pbsv_tbi

        File? sniffles_vcf = CallVariants.sniffles_vcf
        File? sniffles_tbi = CallVariants.sniffles_tbi
        File? sniffles_snf = CallVariants.sniffles_snf

        File? sniffles_phased_vcf = CallVariants.sniffles_phased_vcf
        File? sniffles_phased_tbi = CallVariants.sniffles_phased_tbi
        File? sniffles_phased_snf = CallVariants.sniffles_phased_snf

        File? clair_vcf = CallVariants.clair_vcf
        File? clair_tbi = CallVariants.clair_tbi
        File? clair_gvcf = CallVariants.clair_gvcf
        File? clair_gtbi = CallVariants.clair_gtbi

        File? dv_g_vcf = CallVariants.dv_g_vcf
        File? dv_g_tbi = CallVariants.dv_g_tbi
        File? dv_margin_phased_vcf = CallVariants.dv_margin_phased_vcf
        File? dv_margin_phased_tbi = CallVariants.dv_margin_phased_tbi
        File? dv_vcf_margin_phasing_stats_tsv = CallVariants.dv_vcf_margin_phasing_stats_tsv
        File? dv_vcf_margin_phasing_stats_gtf = CallVariants.dv_vcf_margin_phasing_stats_gtf
        File? dv_whatshap_phased_vcf = CallVariants.dv_whatshap_phased_vcf
        File? dv_whatshap_phased_tbi = CallVariants.dv_whatshap_phased_tbi
        File? dv_vcf_whatshap_phasing_stats_tsv = CallVariants.dv_vcf_whatshap_phasing_stats_tsv
        File? dv_vcf_whatshap_phasing_stats_gtf = CallVariants.dv_vcf_whatshap_phasing_stats_gtf
        String? dv_nongpu_resources_usage_visual = CallVariants.dv_nongpu_resources_usage_visual
    }
}
