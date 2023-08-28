version 1.0

import "../../pipelines/TechAgnostic/Utility/ShardWholeGenome.wdl"  # this isn't optimal; the choice was made assuming ShardWholeGenome could be useful for other users as well.
import "../../deprecated/tasks/PEPPER-MARGIN-DeepVariant.wdl" as PMDV

import "DeepVariant.wdl"
import "Clair.wdl" as Clair3

import "PhaseSmallVariantsAndTagBam.wdl" as PhaseAndTag

workflow Work {
    meta {
        description: "Call small variants using reads-based methods (i.e. not for assembly-contig-based methods)."
    }
    parameter_meta {
        # inputs
        prefix: "Prefix for output files"
        per_chr_bam_bai_and_id: "WGS bam sharded per chromosome/contig."
        is_ont: "If the input data is ONT"
        is_r10_4_pore_or_later: "If the ONT input data is generated on R10.4 simples/duplex pores."
        model_for_dv_andor_pepper: "Model string to be used on DV or the PEPPER-Margin-DeepVariant toolchain. Please refer to their github pages for accepted values."
        ref_scatter_interval_list_locator: "A file holding paths to interval_list files, used for custom sharding the of the input BAM; when not provided, will shard WG by contig (possibly slower)"
        ref_scatter_interval_list_ids: "A file that gives short IDs to the interval_list files; when not provided, will shard WG by contig (possibly slower)"
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
        # sample info
        File bam
        File bai
        String prefix

        Array[Pair[String, Pair[File, File]]] per_chr_bam_bai_and_id

        Boolean is_ont
        Boolean is_r10_4_pore_or_later
        String model_for_dv_andor_pepper

        # reference info
        Map[String, String] ref_map

        File? ref_scatter_interval_list_locator
        File? ref_scatter_interval_list_ids

        # smallVar-specific args
        Boolean run_clair3
        Boolean use_margin_for_tagging

        # optimization
        Int dv_threads
        Int dv_memory
        Boolean use_gpu = false
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    output {
        File? clair_vcf = RunClair3.clair_vcf
        File? clair_tbi = RunClair3.clair_tbi
        File? clair_gvcf = RunClair3.clair_gvcf
        File? clair_gtbi = RunClair3.clair_gtbi

        File haplotagged_bam = use_this_haptag_bam
        File haplotagged_bai = use_this_haptag_bai
        String haplotagged_bam_tagger = use_this_haptagger

        # this block available only for legacy ONT data (those older than R10.4)
        File? legacy_g_vcf              = WorkOnLegacyONTdata.legacy_ont_dvp_g_vcf
        File? legacy_g_tbi              = WorkOnLegacyONTdata.legacy_ont_dvp_g_tbi
        File? legacy_phased_vcf         = WorkOnLegacyONTdata.legacy_ont_dvp_phased_vcf
        File? legacy_phased_tbi         = WorkOnLegacyONTdata.legacy_ont_dvp_phased_tbi
        File? legacy_phasing_stats_tsv  = WorkOnLegacyONTdata.legacy_ont_dvp_phased_vcf_stats_tsv
        File? legacy_phasing_stats_gtf  = WorkOnLegacyONTdata.legacy_ont_dvp_phased_vcf_stats_gtf

        # this block available for CCS and modern ONT data
        File? dv_g_vcf = DV.g_vcf
        File? dv_g_tbi = DV.g_tbi

        File? dv_margin_phased_vcf = PnT.margin_phased_vcf
        File? dv_margin_phased_tbi = PnT.margin_phased_tbi

        File? dv_vcf_margin_phasing_stats_tsv = PnT.margin_phasing_stats_tsv
        File? dv_vcf_margin_phasing_stats_gtf = PnT.margin_phasing_stats_gtf

        File? dv_whatshap_phased_vcf = PnT.whatshap_phased_vcf
        File? dv_whatshap_phased_tbi = PnT.whatshap_phased_tbi

        File? dv_vcf_whatshap_phasing_stats_tsv = PnT.whatshap_phasing_stats_tsv
        File? dv_vcf_whatshap_phasing_stats_gtf = PnT.whatshap_phasing_stats_gtf

        Array[File]? dv_nongpu_resources_usage_log = DV.nongpu_resource_usage_logs
        Array[File]? dv_nongpu_resources_usage_visual = DV.nongpu_resource_usage_visual
    }

    ####################################################################################################################################
    # custom-shard WG (for load-balancing)
    ####################################################################################################################################
    # but if custom sharding isn't requested, then per-chr sharding is already done, so no need to redo
    if (defined(ref_scatter_interval_list_locator)) {
        call ShardWholeGenome.Split as CustomSplitBamForSmallVar {
            input:
                ref_dict = ref_map['dict'],
                bam = bam,
                bai = bai,
                ref_scatter_interval_list_locator = ref_scatter_interval_list_locator,
                ref_scatter_interval_list_ids = ref_scatter_interval_list_ids
        }
    }
    Array[Pair[String, Pair[File, File]]] how_to_shard_wg_for_calling = select_first([CustomSplitBamForSmallVar.id_bam_bai_of_shards,
                                                                                      per_chr_bam_bai_and_id])

    ####################################################################################################################################
    # DV, major workhorse
    ####################################################################################################################################
    if ((!is_ont) || is_r10_4_pore_or_later) { # pacbio or recent ONT data

        call DeepVariant.Run as DV {
            input:
                how_to_shard_wg_for_calling = how_to_shard_wg_for_calling,
                prefix = prefix,
                model_for_dv_andor_pepper = model_for_dv_andor_pepper,

                ref_map = ref_map,

                dv_threads = dv_threads,
                dv_memory = dv_memory,
                use_gpu = use_gpu,
                zones = zones
        }

        call PhaseAndTag.Run as PnT {
            input:
                use_margin_for_tagging = use_margin_for_tagging,

                bam = bam,
                bai = bai,
                per_chr_bam_bai_and_id = per_chr_bam_bai_and_id,
                is_ont = is_ont,

                unphased_vcf = DV.vcf,
                unphased_tbi = DV.tbi,

                ref_map = ref_map,

                zones = zones
        }
    }

    if (is_ont && (!is_r10_4_pore_or_later)) { # legacy (<R10.4) ONT data
        call PMDV.Run as WorkOnLegacyONTdata {
            input:
                how_to_shard_wg_for_calling = how_to_shard_wg_for_calling,
                prefix = prefix,
                model_for_pepper_margin_dv = model_for_dv_andor_pepper,
                ref_map = ref_map,
                dv_threads = dv_threads,
                dv_memory = dv_memory,
                zones = zones
        }
    }

    File use_this_haptag_bam = select_first([PnT.hap_tagged_bam, WorkOnLegacyONTdata.legacy_ont_dvp_haplotagged_bam])
    File use_this_haptag_bai = select_first([PnT.hap_tagged_bai, WorkOnLegacyONTdata.legacy_ont_dvp_haplotagged_bai])
    String use_this_haptagger = select_first([PnT.haplotagged_bam_tagger, "Legacy"])

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
    }
}
