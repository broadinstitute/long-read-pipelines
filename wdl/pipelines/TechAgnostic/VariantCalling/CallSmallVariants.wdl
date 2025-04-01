version 1.0

import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../../structs/ReferenceMetadata.wdl"
import "../../../tasks/Utility/Utils.wdl"

import "../Utility/ShardWholeGenome.wdl"

import "../../../tasks/VariantCalling/DeepVariant.wdl"
import "../../../tasks/VariantCalling/Clair.wdl" as Clair3

import "../Annotation/PhaseSmallVariantsAndTagBam.wdl" as PhaseAndTag

import "../../../deprecated/tasks/PEPPER-MARGIN-DeepVariant.wdl" as PMDV  # this isn't optimal; but sometimes we have legacy data to handle.

struct SmallVarJobConfig {
    String? haploid_contigs

    # optimization
    Int dv_threads
    Int dv_memory
    Boolean use_gpu

    Boolean run_clair3
    Boolean phase_and_tag
    Boolean use_margin_for_tagging

    String? gcp_zones
}

workflow Work {
    meta {
        description: "Call small variants using reads-based methods (i.e. not for assembly-contig-based methods)."
    }
    parameter_meta {
        # inputs
        sex: "biological sex of the sample; accepted value are [F, M, NA]"
        prefix: "Prefix for output files"
        per_chr_bam_bai_and_id: "WGS bam sharded per chromosome/contig."
        is_ont: "If the input data is ONT"
        is_r10_4_pore_or_later: "If the ONT input data is generated on R10.4 simples/duplex pores."
        model_for_dv_andor_pepper: "Model string to be used on DV or the PEPPER-Margin-DeepVariant toolchain. Please refer to their github pages for accepted values."
        phase_and_tag: "if turned on, small variants will be phased using both WhatsHap and Margin, then the BAM will be haplotagged with either the output of WhatsHap or Margin (depending on use_margin_for_tagging); then the haplotagged BAM will be used for calling phased-SV again with Sniffles2. Obviously, this prolongs the runtime significantly. For ONT data on pores older than R10.4, having this turned off means phased VCF and haplotagged BAM will not be output."
        haploid_contigs: "Optimization since DV 1.6 to improve calling on haploid contigs (e.g. human allosomes); see DV github page for more info."
        use_gpu: "Use GPU acceleration for DV (or PEPPER) or not"
        use_margin_for_tagging: "if false, will use margin-phased VCF for haplotagging the BAM; applicable only when input data isn't ONT data with pore older than R10.4"

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

        dv_nongpu_resources_usage_log: "Resource usage monitoring log for DV (per shard); available for CCS data and ONT data generated with pores >= R10.4."
        dv_nongpu_resources_usage_visual: "Resource usage monitoring log visualization for DV (per shard); available for CCS data and ONT data generated with pores >= R10.4."
    }
    input {
        String gcs_variants_out_dir
        String gcs_tagged_bam_out_dir

        # sample info
        File bam
        File bai
        String sex
        String prefix

        Array[Pair[String, Pair[File, File]]] per_chr_bam_bai_and_id
        Boolean force_per_chr_sharding_scheme = false

        Boolean is_ont
        Boolean is_r10_4_pore_or_later
        String model_for_dv_andor_pepper

        # reference info
        File ref_bundle_json_file

        # phasing and read-haplotaging desired or not
        Boolean phase_and_tag

        # smallVar-specific args
        Boolean run_clair3
        Boolean use_margin_for_tagging
        String? haploid_contigs

        # optimization
        Int dv_threads
        Int dv_memory
        Boolean use_gpu = false
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        # DO NOT EVER DEFINE THESE TWO VARIABLES; THEY ARE USED AS A WDL TRICK. THEY SHOULD NEVER BE DEFINED.
        String? dummy_a
        File? dummy_b
    }
    output {
        File? clair_vcf  = FinalizeClairVcf.gcs_path
        File? clair_tbi  = FinalizeClairTbi.gcs_path
        File? clair_gvcf = FinalizeClairGVcf.gcs_path
        File? clair_gtbi = FinalizeClairGTbi.gcs_path

        File? haplotagged_bam = use_this_haptag_bam
        File? haplotagged_bai = use_this_haptag_bai
        String? haplotagged_bam_tagger = use_this_haptagger

        # this block available only for legacy ONT data (those older than R10.4)
        File? legacy_g_vcf              = FinalizeLegacyGVcf.gcs_path
        File? legacy_g_tbi              = FinalizeLegacyGTbi.gcs_path
        File? legacy_phased_vcf         = FinalizeLegacyPhasedVcf.gcs_path
        File? legacy_phased_tbi         = FinalizeLegacyPhasedTbi.gcs_path
        File? legacy_phasing_stats_tsv  = FinalizeLegacyPhaseStatsTSV.gcs_path
        File? legacy_phasing_stats_gtf  = FinalizeLegacyPhaseStatsGTF.gcs_path

        # this block available for CCS and modern ONT data
        File? dv_g_vcf = FinalizeDVgVcf.gcs_path
        File? dv_g_tbi = FinalizeDVgTbi.gcs_path

        String? dv_native_visual_report_html = FinalizeDVNativeVisualReportHTMLs.gcs_dir
        String? dv_nongpu_resources_usage_log = FinalizeDVResourceUsagesLogging.gcs_dir
        String? dv_nongpu_resources_usage_visual = FinalizeDVResourceUsagesVisual.gcs_dir

        File? dv_margin_phased_vcf = PhaseThenTag.margin_phased_vcf
        File? dv_margin_phased_tbi = PhaseThenTag.margin_phased_tbi

        File? dv_vcf_margin_phasing_stats_tsv = PhaseThenTag.margin_phasing_stats_tsv
        File? dv_vcf_margin_phasing_stats_gtf = PhaseThenTag.margin_phasing_stats_gtf

        File? dv_whatshap_phased_vcf = PhaseThenTag.whatshap_phased_vcf
        File? dv_whatshap_phased_tbi = PhaseThenTag.whatshap_phased_tbi

        File? dv_vcf_whatshap_phasing_stats_tsv = PhaseThenTag.whatshap_phasing_stats_tsv
        File? dv_vcf_whatshap_phasing_stats_gtf = PhaseThenTag.whatshap_phasing_stats_gtf
    }
    if (defined(dummy_a) || defined(dummy_b)) {
        call Utils.StopWorkflow as DummyVariablesAreAccidentallyDefined { input:
            reason = "dummy_a and dummy_b are used as a WDL trick; they should never be defined."
        }
    }

    HumanReferenceBundle ref_bundle = read_json(ref_bundle_json_file)
    Map[String, String] ref_map = ref_bundle
    ####################################################################################################################################
    # custom-shard WG (for load-balancing)
    ####################################################################################################################################
    # but if custom sharding isn't requested, then per-chr sharding is already done, so no need to redo
    if (defined(ref_bundle.size_balanced_scatter_intervallists_locators) && (!force_per_chr_sharding_scheme)) {
        call ShardWholeGenome.Split as CustomSplitBamForSmallVar {
            input:
                bam = bam,
                bai = bai,
                ref_dict = ref_bundle.dict,
                ref_scatter_interval_list_ids = select_first([ref_bundle.size_balanced_scatter_interval_ids]),
                ref_scatter_interval_list_locator = select_first([ref_bundle.size_balanced_scatter_intervallists_locators])
        }
    }
    Array[Pair[String, Pair[File, File]]] how_to_shard_wg_for_calling = if (force_per_chr_sharding_scheme)
                                                                        then per_chr_bam_bai_and_id
                                                                        else select_first([CustomSplitBamForSmallVar.id_bam_bai_of_shards, per_chr_bam_bai_and_id])

    ####################################################################################################################################
    # DV, major workhorse
    ####################################################################################################################################
    Boolean is_legacy_ont = is_ont && (!is_r10_4_pore_or_later)

    if (!is_legacy_ont) { # pacbio or recent ONT data
        call DeepVariant.Run as DV {
            input:
                how_to_shard_wg_for_calling = how_to_shard_wg_for_calling,

                prefix = prefix,
                model_for_dv_andor_pepper = model_for_dv_andor_pepper,

                ref_bundle = ref_bundle,

                haploid_contigs = if (sex == "M") then haploid_contigs    else dummy_a,
                par_regions_bed = if (sex == "M") then ref_bundle.PAR_bed else dummy_b,

                dv_threads = dv_threads,
                dv_memory = dv_memory,
                use_gpu = use_gpu,
                zones = zones
        }
        call FF.FinalizeToFile as FinalizeDVgVcf { input: outdir = gcs_variants_out_dir, file = DV.g_vcf }
        call FF.FinalizeToFile as FinalizeDVgTbi { input: outdir = gcs_variants_out_dir, file = DV.g_tbi }
        call FF.FinalizeToDir as FinalizeDVNativeVisualReportHTMLs {
            input: files = DV.native_visual_report_htmls, outdir = gcs_variants_out_dir + "/DV_monitoring"
        }
        call FF.FinalizeToDir as FinalizeDVResourceUsagesLogging {
            input: files = DV.nongpu_resource_usage_logs, outdir = gcs_variants_out_dir + "/DV_monitoring"
        }
        call FF.FinalizeToDir as FinalizeDVResourceUsagesVisual {
            input: files = DV.nongpu_resource_usage_visual, outdir = gcs_variants_out_dir + "/DV_monitoring"
        }

        if (phase_and_tag) {
            call PhaseAndTag.Run as PhaseThenTag {
                input:
                    use_margin_for_tagging = use_margin_for_tagging,

                    bam = bam,
                    bai = bai,
                    per_chr_bam_bai_and_id = per_chr_bam_bai_and_id,
                    is_ont = is_ont,

                    unphased_vcf = DV.vcf,
                    unphased_tbi = DV.tbi,

                    ref_map = ref_map,

                    zones = zones,
                    gcs_variants_out_dir = gcs_variants_out_dir,
                    gcs_tagged_bam_out_dir = gcs_tagged_bam_out_dir
            }
        }
    }

    if (is_legacy_ont) { # legacy (<R10.4) ONT data
        call PMDV.Run as WorkOnLegacyONTdata {
            input:
                how_to_shard_wg_for_calling = how_to_shard_wg_for_calling,
                prefix = prefix,
                model_for_pepper_margin_dv = model_for_dv_andor_pepper,
                ref_map = ref_map,
                dv_threads = dv_threads,
                dv_memory = dv_memory,
                phase_and_tag = phase_and_tag,
                zones = zones
        }
        call FF.FinalizeToFile as FinalizeLegacyGVcf { input: outdir = gcs_variants_out_dir, file = WorkOnLegacyONTdata.legacy_ont_dvp_g_vcf }
        call FF.FinalizeToFile as FinalizeLegacyGTbi { input: outdir = gcs_variants_out_dir, file = WorkOnLegacyONTdata.legacy_ont_dvp_g_tbi }
        if (phase_and_tag) {
            call FF.FinalizeToFile as FinalizeLegacyPhasedVcf       { input: outdir = gcs_variants_out_dir, file = select_first([WorkOnLegacyONTdata.legacy_ont_dvp_phased_vcf]) }
            call FF.FinalizeToFile as FinalizeLegacyPhasedTbi       { input: outdir = gcs_variants_out_dir, file = select_first([WorkOnLegacyONTdata.legacy_ont_dvp_phased_tbi]) }
            call FF.FinalizeToFile as FinalizeLegacyPhaseStatsTSV   { input: outdir = gcs_variants_out_dir, file = select_first([WorkOnLegacyONTdata.legacy_ont_dvp_phased_vcf_stats_tsv]) }
            call FF.FinalizeToFile as FinalizeLegacyPhaseStatsGTF   { input: outdir = gcs_variants_out_dir, file = select_first([WorkOnLegacyONTdata.legacy_ont_dvp_phased_vcf_stats_gtf]) }

            call FF.FinalizeToFile as FinalizeLegacyHapTaggedBam { input: outdir = gcs_tagged_bam_out_dir, file = select_first([WorkOnLegacyONTdata.legacy_ont_dvp_haplotagged_bam]) }
            call FF.FinalizeToFile as FinalizeLegacyHapTaggedBai { input: outdir = gcs_tagged_bam_out_dir, file = select_first([WorkOnLegacyONTdata.legacy_ont_dvp_haplotagged_bai]) }
        }
    }

    if (phase_and_tag) {
        File use_this_haptag_bam  = select_first([PhaseThenTag.hap_tagged_bam, FinalizeLegacyHapTaggedBam.gcs_path])
        File use_this_haptag_bai  = select_first([PhaseThenTag.hap_tagged_bai, FinalizeLegacyHapTaggedBai.gcs_path])
        String use_this_haptagger = select_first([PhaseThenTag.haplotagged_bam_tagger, "Legacy"])
    }

    ####################################################################################################################################
    # run clair3, if so requested
    ####################################################################################################################################
    if (run_clair3) {
        call Clair3.Run as RunClair3 {
            input:
                how_to_shard_wg_for_calling = how_to_shard_wg_for_calling,
                is_ont = is_ont, prefix = prefix,
                ref_map = ref_map, zones = zones
        }

        call FF.FinalizeToFile as FinalizeClairVcf  { input: outdir = gcs_variants_out_dir, file = RunClair3.clair_vcf }
        call FF.FinalizeToFile as FinalizeClairTbi  { input: outdir = gcs_variants_out_dir, file = RunClair3.clair_tbi }
        call FF.FinalizeToFile as FinalizeClairGVcf { input: outdir = gcs_variants_out_dir, file = RunClair3.clair_gvcf }
        call FF.FinalizeToFile as FinalizeClairGTbi { input: outdir = gcs_variants_out_dir, file = RunClair3.clair_gtbi }
    }
}
