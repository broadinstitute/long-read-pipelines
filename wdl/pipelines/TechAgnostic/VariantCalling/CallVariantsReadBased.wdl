version 1.0

import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../../tasks/Utility/Utils.wdl"

import "../../../structs/ReferenceMetadata.wdl"

import "../Utility/ShardWholeGenome.wdl"

import "CallStructuralVariants.wdl"
import "CallSmallVariants.wdl"

import "../../../tasks/VariantCalling/Sniffles2.wdl"


workflow CallVariants {

    meta {
        description: "A workflow for calling small and/or structural variants from an aligned BAM file. Note this calls out to read-based methods, not assembly-based methods. This also does not support CLR data."
    }

    parameter_meta {
        bam: "Aligned BAM file"
        bai: "Index for the aligned BAM file"
        sex: "biological sex of the sample; accepted value are [F, M, NA]"
        prefix: "Prefix for output files"

        is_ont: "If the input data is generated on the ONT platform"
        is_r10_4_pore_or_later: "tell us which pore version was used to generate the data. When true, will use the DV (>=1.5.0) toolchain."
        model_for_dv_andor_pepper: "model string to be used on DV or the PEPPER-Margin-DeepVariant toolchain. Please refer to their github pages for accepted values."

        ref_bundle_json_file: "a json file holding reference file location and auxillary file locations; see HumanReferenceBundle struct defined in ReferenceMetadata"

        small_variant_calling_options_json: "a json file holding config for small variant calling (see struct SmallVarJobConfig in CallSmallVariants.wdl for detail; when omitted, will skip small variant calling"
        sv_calling_options_json: "a json file holding config for SV calling (see struct SVCallingConfig in CallStructuralVariants.wdl for detail; when omitted, will skip SV calling"

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
        String gcs_out_dir

        # sample info
        File bam
        File bai
        String sex
        String prefix

        # data type info
        Boolean is_ont
        Boolean is_r10_4_pore_or_later
        String model_for_dv_andor_pepper

        # reference-specific
        File ref_bundle_json_file

        File? small_variant_calling_options_json
        File? sv_calling_options_json

        Boolean force_per_chr_sharding_scheme = false

        Array[String] gcp_zones = ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]
    }

    if ((!defined(sv_calling_options_json)) && (!defined(small_variant_calling_options_json))) {
        call Utils.StopWorkflow { input: reason = "Why are you calling me if your want neither small variants nor SVs?"}
    }

    ######################################################################
    # Block for prepping inputs
    ######################################################################
    HumanReferenceBundle ref_bundle = read_json(ref_bundle_json_file)

    call GU.CollapseArrayOfStrings as get_zones {input: input_array = gcp_zones, joiner = " "}
    String wdl_parsable_zones = get_zones.collapsed

    # needed for whatshap phasing anyway, so this can be used by SV calling
    call ShardWholeGenome.Split as SplitBamByChr { input: ref_dict = ref_bundle.dict, bam = bam, bai = bai, }

    ######################################################################
    # Block for small variants handling
    ######################################################################
    if (defined(small_variant_calling_options_json)) {

        SmallVarJobConfig snp_options = read_json(select_first([small_variant_calling_options_json]))

        call CallSmallVariants.Work as SmallVarJob {
            input:
                bam = bam,
                bai = bai,
                sex = sex,
                prefix = prefix,

                per_chr_bam_bai_and_id = SplitBamByChr.id_bam_bai_of_shards,
                force_per_chr_sharding_scheme = force_per_chr_sharding_scheme,

                is_ont = is_ont,
                is_r10_4_pore_or_later = is_r10_4_pore_or_later,
                model_for_dv_andor_pepper = model_for_dv_andor_pepper,

                gcs_variants_out_dir = sub(gcs_out_dir, "/$", "") + "/variants/small",
                gcs_tagged_bam_out_dir = sub(gcs_out_dir, "/$", "") + "/alignments",

                ref_bundle_json_file = ref_bundle_json_file,

                run_clair3 = snp_options.run_clair3,

                phase_and_tag = snp_options.phase_and_tag,
                use_margin_for_tagging = snp_options.use_margin_for_tagging,

                haploid_contigs = snp_options.haploid_contigs,

                dv_threads = snp_options.dv_threads,
                dv_memory = snp_options.dv_memory,
                use_gpu = snp_options.use_gpu,
                zones = wdl_parsable_zones
        }
    }

    ######################################################################
    # Block for SV handling
    ######################################################################
    String sv_dir = sub(gcs_out_dir, "/$", "") + "/variants/sv"
    if (defined(sv_calling_options_json)) {

        SVCallingConfig sv_options = read_json(select_first([sv_calling_options_json]))

        call CallStructuralVariants.Work as SVjob {
            input:
                is_hifi = !is_ont,
                is_ont = is_ont,

                bam = bam,
                bai = bai,
                prefix = prefix,

                per_chr_bam_bai_and_id = SplitBamByChr.id_bam_bai_of_shards,

                gcs_out_dir = sv_dir,

                ref_bundle_json_file = ref_bundle_json_file,

                minsvlen = sv_options.min_sv_len,
                pbsv_discover_per_chr = sv_options.pbsv_discover_per_chr,

                zones = wdl_parsable_zones
        }
    }

    ######################################################################
    # Experiment with Sniffles-2 phased SV calling
    ######################################################################
    Boolean call_both = defined(sv_calling_options_json) && defined(small_variant_calling_options_json)
    if (call_both) {  # but do so lazily
        SmallVarJobConfig opt_a = read_json(select_first([small_variant_calling_options_json]))
        SVCallingConfig opt_b = read_json(select_first([sv_calling_options_json]))
        if (opt_a.phase_and_tag) {
            File m = select_first([SmallVarJob.haplotagged_bam])
            File i = select_first([SmallVarJob.haplotagged_bai])
            call Utils.InferSampleName { input: bam = m, bai = i }
            call Sniffles2.SampleSV as SnifflesPhaseSV {
                input:
                    bam = m, bai = i, sample_id = InferSampleName.sample_name,
                    prefix = prefix, tandem_repeat_bed = ref_bundle.tandem_repeat_bed,
                    minsvlen = opt_b.min_sv_len,
                    phase_sv = true
            }
            call FF.FinalizeToFile as FinalizePhasedSnifflesVcf { input: outdir = sv_dir, file = SnifflesPhaseSV.vcf }
            call FF.FinalizeToFile as FinalizePhasedSnifflesTbi { input: outdir = sv_dir, file = SnifflesPhaseSV.tbi }
            call FF.FinalizeToFile as FinalizePhasedSnifflesSnf { input: outdir = sv_dir, file = SnifflesPhaseSV.snf }
        }
    }

    output {
        ####### SVs
        File? pbsv_vcf = SVjob.pbsv_vcf
        File? pbsv_tbi = SVjob.pbsv_tbi

        File? sniffles_vcf = SVjob.sniffles_vcf
        File? sniffles_tbi = SVjob.sniffles_tbi
        File? sniffles_snf = SVjob.sniffles_snf

        File? sniffles_phased_vcf = FinalizePhasedSnifflesVcf.gcs_path
        File? sniffles_phased_tbi = FinalizePhasedSnifflesTbi.gcs_path
        File? sniffles_phased_snf = FinalizePhasedSnifflesSnf.gcs_path

        ####### small vars
        # clair
        File? clair_vcf = SmallVarJob.clair_vcf
        File? clair_tbi = SmallVarJob.clair_tbi
        File? clair_gvcf = SmallVarJob.clair_gvcf
        File? clair_gtbi = SmallVarJob.clair_gtbi

        # tagging
        File? haplotagged_bam = SmallVarJob.haplotagged_bam
        File? haplotagged_bai = SmallVarJob.haplotagged_bai
        String? haplotagged_bam_tagger = SmallVarJob.haplotagged_bam_tagger

        # available for CCS and ONT >= R10.4 data, if small variants are requested
        File? dv_g_vcf = SmallVarJob.dv_g_vcf
        File? dv_g_tbi = SmallVarJob.dv_g_tbi
        File? dv_margin_phased_vcf = SmallVarJob.dv_margin_phased_vcf
        File? dv_margin_phased_tbi = SmallVarJob.dv_margin_phased_tbi
        File? dv_vcf_margin_phasing_stats_tsv = SmallVarJob.dv_vcf_margin_phasing_stats_tsv
        File? dv_vcf_margin_phasing_stats_gtf = SmallVarJob.dv_vcf_margin_phasing_stats_gtf
        File? dv_whatshap_phased_vcf = SmallVarJob.dv_whatshap_phased_vcf
        File? dv_whatshap_phased_tbi = SmallVarJob.dv_whatshap_phased_tbi
        File? dv_vcf_whatshap_phasing_stats_tsv = SmallVarJob.dv_vcf_whatshap_phasing_stats_tsv
        File? dv_vcf_whatshap_phasing_stats_gtf = SmallVarJob.dv_vcf_whatshap_phasing_stats_gtf
        String? dv_nongpu_resources_usage_visual = SmallVarJob.dv_nongpu_resources_usage_visual
        String? dv_native_visual_report_html = SmallVarJob.dv_native_visual_report_html

        # available for ONT < R10.4 data, if small variants are requested
        File? legacy_g_vcf = SmallVarJob.legacy_g_vcf
        File? legacy_g_tbi = SmallVarJob.legacy_g_tbi
        File? legacy_phased_vcf = SmallVarJob.legacy_phased_vcf
        File? legacy_phased_tbi = SmallVarJob.legacy_phased_tbi
        File? legacy_phasing_stats_tsv = SmallVarJob.legacy_phasing_stats_tsv
        File? legacy_phasing_stats_gtf = SmallVarJob.legacy_phasing_stats_gtf
    }
}
