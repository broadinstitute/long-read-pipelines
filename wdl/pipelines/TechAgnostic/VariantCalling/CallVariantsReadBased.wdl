version 1.0

import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../../tasks/Utility/Utils.wdl"

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
        prefix: "Prefix for output files"

        is_ont: "If the input data is generated on the ONT platform"
        is_r10_4_pore_or_later: "tell us which pore version was used to generate the data. When true, will use the DV (>=1.5.0) toolchain."
        model_for_dv_andor_pepper: "model string to be used on DV or the PEPPER-Margin-DeepVariant toolchain. Please refer to their github pages for accepted values."

        ref_map_file: "table indicating reference sequence and auxillary file locations"
        ref_scatter_interval_list_locator: "A file holding paths to interval_list files, used for custom sharding the of the input BAM; when not provided, will shard WG by contig (possibly slower)"
        ref_scatter_interval_list_ids: "A file that gives short IDs to the interval_list files; when not provided, will shard WG by contig (possibly slower)"

        phase_and_tag: "if turned on, small variants will be phased using both WhatsHap and Margin, then the BAM will be haplotagged with either the output of WhatsHap or Margin (depending on use_margin_for_tagging); then the haplotagged BAM will be used for calling phased-SV again with Sniffles2. Obviously, this prolongs the runtime significantly. Has no effect on ONT data on pores older than R10.4."

        call_svs: "Call structural variants or not"
        minsvlen: "Minimum SV length in bp; only affects Sniffles 2 calls."
        pbsv_discover_per_chr: "Run the discover stage of PBSV per chromosome"

        call_small_variants: "Call small variants or not"
        run_clair3: "to turn on Clair3 analysis or not (non-trivial increase in cost and runtime)"
        use_margin_for_tagging: "if false, will use margin-phased small-variant VCF for haplotagging the BAM; applicable only when input data isn't ONT data with pore older than R10.4"

        dv_threads: "number of threads for DeepVariant"
        dv_memory:  "memory for DeepVariant"
        use_gpu: "to use GPU acceleration or not on DeepVariant"

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
        String prefix

        # data type info
        Boolean is_ont
        Boolean is_r10_4_pore_or_later
        String model_for_dv_andor_pepper

        # reference-specific
        File ref_map_file
        File? ref_scatter_interval_list_locator
        File? ref_scatter_interval_list_ids

        # sv-specific args
        Boolean call_svs
        Boolean pbsv_discover_per_chr
        Int minsvlen = 20

        # smallVar-specific args
        Boolean call_small_variants
        Boolean run_clair3
        Boolean use_margin_for_tagging

        # phasing and read-haplotaging desired or not
        Boolean phase_and_tag

        # optimization, balancing between throughput, wallclock time, and cost
        Int dv_threads
        Int dv_memory
        Boolean use_gpu = false
        Array[String] gcp_zones = ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]
    }

    if ((!call_svs) && (!call_small_variants)) {
        call Utils.StopWorkflow { input: reason = "Why are you calling me if your want neither small variants nor SVs?"}
    }

    ######################################################################
    # Block for prepping inputs
    ######################################################################
    Map[String, String] ref_map = read_map(ref_map_file)

    call GU.CollapseArrayOfStrings as get_zones {input: input_array = gcp_zones, joiner = " "}
    String wdl_parsable_zones = get_zones.collapsed

    # needed for whatshap phasing anyway, so this can be used by SV calling
    call ShardWholeGenome.Split as SplitBamByChr { input: ref_dict = ref_map['dict'], bam = bam, bai = bai, }

    ######################################################################
    # Block for small variants handling
    ######################################################################
    if (call_small_variants) {

        call CallSmallVariants.Work as SmallVarJob {
            input:
                bam = bam,
                bai = bai,
                prefix = prefix,

                per_chr_bam_bai_and_id = SplitBamByChr.id_bam_bai_of_shards,

                is_ont = is_ont,
                is_r10_4_pore_or_later = is_r10_4_pore_or_later,
                model_for_dv_andor_pepper = model_for_dv_andor_pepper,

                ref_map = ref_map,
                ref_scatter_interval_list_locator = ref_scatter_interval_list_locator,
                ref_scatter_interval_list_ids = ref_scatter_interval_list_ids,

                run_clair3 = run_clair3,

                phase_and_tag = phase_and_tag,
                use_margin_for_tagging = use_margin_for_tagging,

                dv_threads = dv_threads,
                dv_memory = dv_memory,
                use_gpu = use_gpu,
                zones = wdl_parsable_zones,

                gcs_variants_out_dir = sub(gcs_out_dir, "/$", "") + "/variants/small",
                gcs_tagged_bam_out_dir = sub(gcs_out_dir, "/$", "") + "/alignments"
        }
    }

    ######################################################################
    # Block for SV handling
    ######################################################################
    String sv_dir = sub(gcs_out_dir, "/$", "") + "/variants/sv"
    if (call_svs) {
        call CallStructuralVariants.Work as SVjob {
            input:
                is_hifi = !is_ont,
                is_ont = is_ont,

                bam = bam,
                bai = bai,
                prefix = prefix,

                per_chr_bam_bai_and_id = SplitBamByChr.id_bam_bai_of_shards,

                ref_map = ref_map,

                minsvlen = minsvlen,

                pbsv_discover_per_chr = pbsv_discover_per_chr,

                zones = wdl_parsable_zones,
                gcs_out_dir = sv_dir
        }
    }

    ######################################################################
    # Experiment with Sniffles-2 phased SV calling
    ######################################################################
    if (call_svs && call_small_variants && phase_and_tag) {  # but do so lazily
        File m = select_first([SmallVarJob.haplotagged_bam])
        File i = select_first([SmallVarJob.haplotagged_bai])
        call Utils.InferSampleName { input: bam = m, bai = i }
        call Sniffles2.SampleSV as SnifflesPhaseSV {
            input:
                bam = m, bai = i, sample_id = InferSampleName.sample_name,
                prefix = prefix, tandem_repeat_bed = ref_map['tandem_repeat_bed'],
                minsvlen = minsvlen,
                phase_sv = true
        }
        call FF.FinalizeToFile as FinalizePhasedSnifflesVcf { input: outdir = sv_dir, file = SnifflesPhaseSV.vcf }
        call FF.FinalizeToFile as FinalizePhasedSnifflesTbi { input: outdir = sv_dir, file = SnifflesPhaseSV.tbi }
        call FF.FinalizeToFile as FinalizePhasedSnifflesSnf { input: outdir = sv_dir, file = SnifflesPhaseSV.snf }
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

        # available for ONT < R10.4 data, if small variants are requested
        File? legacy_g_vcf = SmallVarJob.legacy_g_vcf
        File? legacy_g_tbi = SmallVarJob.legacy_g_tbi
        File? legacy_phased_vcf = SmallVarJob.legacy_phased_vcf
        File? legacy_phased_tbi = SmallVarJob.legacy_phased_tbi
        File? legacy_phasing_stats_tsv = SmallVarJob.legacy_phasing_stats_tsv
        File? legacy_phasing_stats_gtf = SmallVarJob.legacy_phasing_stats_gtf
    }
}
