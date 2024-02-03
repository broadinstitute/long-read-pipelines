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

        ref_bundle_json_file: "a json file holding reference file location and auxillary file locations; see HumanReferenceBundle struct defined in ReferenceMetadata"

        bed_to_compute_coverage: "BED file holding regions-of-interest for computing coverage over."
        bed_descriptor: "Description of the BED file, will be used in the file name so be careful naming things"

        small_variant_calling_options_json: "a json file holding config for small variant calling (see struct SmallVarJobConfig in CallSmallVariants.wdl for detail; when omitted, will skip small variant calling"
        sv_calling_options_json: "a json file holding config for SV calling (see struct SVCallingConfig in CallStructuralVariants.wdl for detail; when omitted, will skip SV calling"

        gcp_zones: "which Google Cloud Zone to use (this has implications on how many GPUs are available and egress costs, so configure carefully)"

        # outputs
        haplotagged_bam: "BAM haplotagged using a small variant single-sample VCF."
        haplotagged_bai: "Index for haplotagged_bam."
        haplotagged_bam_tagger: "VCF used for doing the haplotagging. 'Legacy' if the input is ONT data generated on pores before R10.4."

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

        # reference-specific
        File ref_bundle_json_file
        File? bed_to_compute_coverage
        String? bed_descriptor

        # user choice
        File? small_variant_calling_options_json
        File? sv_calling_options_json

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

            is_ont = false,
            bams_suspected_to_contain_dup_record = false,

            bed_to_compute_coverage = bed_to_compute_coverage,
            bed_descriptor = bed_descriptor
    }

    ###########################################################
    if (defined(sv_calling_options_json) || defined(small_variant_calling_options_json)) {
        call VAR.CallVariants {
            input:
                gcs_out_dir = outdir,
                bam    = MergeAndMetrics.aligned_bam,
                bai    = MergeAndMetrics.aligned_bai,
                prefix = sample_name,

                is_ont = false,
                is_r10_4_pore_or_later = false,
                model_for_dv_andor_pepper = 'PACBIO',

                ref_bundle_json_file = ref_bundle_json_file,

                small_variant_calling_options_json = small_variant_calling_options_json,
                sv_calling_options_json = sv_calling_options_json,

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
    }
}
