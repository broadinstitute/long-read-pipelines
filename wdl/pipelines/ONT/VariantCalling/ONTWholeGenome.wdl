version 1.0

import "../../../tasks/Utility/GeneralUtils.wdl" as GU

import "../../TechAgnostic/Utility/MergeSampleBamsAndCollectMetrics.wdl" as MERGE
import "../../TechAgnostic/VariantCalling/CallVariantsReadBased.wdl" as VAR

workflow ONTWholeGenome {

    meta {
        description: "A workflow that performs single sample variant calling on Oxford Nanopore reads from one or more flow cells. The workflow merges multiple flowcells into a single BAM prior to variant calling."
    }
    parameter_meta {
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"

        aligned_bams:       "GCS path to aligned BAM files"
        aligned_bais:       "GCS path to aligned BAM file indices"
        sample_name:        "sample name as encoded in the bams"

        bams_suspected_to_contain_dup_record: "Some ONT output files from basecall dirs have a strange duplicate issue."
        is_r10_4_pore_or_later: "tell us which pore version was used to generate the data. When true, will use the DV (>=1.5.0) toolchain."
        model_for_dv_andor_pepper: "model string to be used on DV or the PEPPER-Margin-DeepVariant toolchain. Please refer to their github pages for accepted values."

        ref_map_file: "table indicating reference sequence and auxillary file locations"

        bed_to_compute_coverage: "BED file holding regions-of-interest for computing coverage over."
        bed_descriptor: "Description of the BED file, will be used in the file name so be careful naming things"

        small_variant_calling_options_json: "a json file holding config for small variant calling (see struct SmallVarJobConfig in CallSmallVariants.wdl for detail; when omitted, will skip small variant calling"
        sv_calling_options_json: "a json file holding config for SV calling (see struct SVCallingConfig in CallStructuralVariants.wdl for detail; when omitted, will skip SV calling"

        gcp_zones: "which Google Cloud Zone to use (this has implications on how many GPUs are available and egress costs, so configure carefully)"

        # outputs
        haplotagged_bam: "BAM haplotagged using a small variant single-sample VCF."
        haplotagged_bai: "Index for haplotagged_bam."
        haplotagged_bam_tagger: "VCF used for doing the haplotagging. 'Legacy' if the input is ONT data generated on pores before R10.4."

        legacy_g_vcf: "PEPPER-MARGIN-DeepVariant gVCF; available only when input is ONT data generated on pores older than R10.4."
        legacy_g_tbi: "Index for PEPPER-MARGIN-DeepVariant gVCF; available only when input is ONT data generated on pores older than R10.4."
        legacy_phased_vcf: "Phased PEPPER-MARGIN-DeepVariant VCF; available only when input is ONT data generated on pores older than R10.4."
        legacy_phased_tbi: "Indes for phased PEPPER-MARGIN-DeepVariant VCF; available only when input is ONT data generated on pores older than R10.4."
        legacy_phasing_stats_tsv: "Phasing stats of legacy_phased_vcf in TSV format; available only when input is ONT data generated on pores older than R10.4."
        legacy_phasing_stats_gtf: "Phasing stats of legacy_phased_vcf in GTF format; available only when input is ONT data generated on pores older than R10.4."

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
        Boolean bams_suspected_to_contain_dup_record
        Boolean is_r10_4_pore_or_later
        String model_for_dv_andor_pepper

        # reference-specific
        File ref_map_file
        File? bed_to_compute_coverage
        String? bed_descriptor

        # user choice
        File? small_variant_calling_options_json
        File? sv_calling_options_json

        Array[String] gcp_zones = ['us-central1-a', 'us-central1-b', 'us-central1-c', 'us-central1-f']
    }

    String workflow_name = "ONTWholeGenome"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_name}/~{sample_name}"

    ###########################################################
    call MERGE.Work as MergeAndMetrics {
        input:
            gcs_out_dir = outdir,

            sample_name = sample_name,
            aligned_bams = aligned_bams,
            aligned_bais = aligned_bais,

            is_ont = true,
            bams_suspected_to_contain_dup_record = bams_suspected_to_contain_dup_record,

            bed_to_compute_coverage = bed_to_compute_coverage,
            bed_descriptor = bed_descriptor
    }

    ####################################################################################################
    if (defined(sv_calling_options_json) || defined(small_variant_calling_options_json)) {
        call VAR.CallVariants {
            input:
                gcs_out_dir = outdir,

                bam    = MergeAndMetrics.aligned_bam,
                bai    = MergeAndMetrics.aligned_bai,
                prefix = sample_name,

                is_ont = true,
                is_r10_4_pore_or_later = is_r10_4_pore_or_later,
                model_for_dv_andor_pepper = model_for_dv_andor_pepper,

                ref_map_file = ref_map_file,

                small_variant_calling_options_json = small_variant_calling_options_json,
                sv_calling_options_json = sv_calling_options_json,

                gcp_zones = gcp_zones
        }
    }

    ####################################################################################################
    call GU.GetTodayDate as today {}

    ###########################################################
    output {
        String last_processing_date = today.yyyy_mm_dd

        ########################################
        File aligned_bam = MergeAndMetrics.aligned_bam
        File aligned_bai = MergeAndMetrics.aligned_bai

        Float coverage = MergeAndMetrics.coverage
        File? bed_cov_summary = MergeAndMetrics.bed_cov_summary

        Map[String, Float] alignment_metrics = MergeAndMetrics.alignment_metrics

        ########################################
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

        # available for ONT >= R10.4 data, if small variants are requested
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
        String? dv_native_visual_report_html = CallVariants.dv_native_visual_report_html

        # available for ONT < R10.4 data, if small variants are requested
        File? legacy_g_vcf = CallVariants.legacy_g_vcf
        File? legacy_g_tbi = CallVariants.legacy_g_tbi
        File? legacy_phased_vcf = CallVariants.legacy_phased_vcf
        File? legacy_phased_tbi = CallVariants.legacy_phased_tbi
        File? legacy_phasing_stats_tsv = CallVariants.legacy_phasing_stats_tsv
        File? legacy_phasing_stats_gtf = CallVariants.legacy_phasing_stats_gtf
    }
}
