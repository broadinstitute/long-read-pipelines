version 1.0

import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../../tasks/Utility/Utils.wdl"

import "../Utility/ShardWholeGenome.wdl"

import "../../../tasks/VariantCalling/CallStructuralVariants.wdl"
import "../../../tasks/VariantCalling/CallSmallVariants.wdl"

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

        call_svs: "Call structural variants or not"
        minsvlen: "Minimum SV length in bp (default: 50)"
        pbsv_discover_per_chr: "Run the discover stage of PBSV per chromosome"

        call_small_variants: "Call small variants or not"
        run_clair3: "to turn on Clair3 analysis or not (non-trivial increase in cost and runtime)"
        dv_threads: "number of threads for DeepVariant"
        dv_memory:  "memory for DeepVariant"
        use_gpu: "to use GPU acceleration or not on DeepVariant"
    }

    input {
        String? gcs_out_dir

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
        Int minsvlen = 50

        # smallVar-specific args
        Boolean call_small_variants
        Boolean run_clair3

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

                ref_fasta = ref_map['fasta'],
                ref_fasta_fai = ref_map['fai'],
                ref_dict = ref_map['dict'],

                ref_scatter_interval_list_locator = ref_scatter_interval_list_locator,
                ref_scatter_interval_list_ids = ref_scatter_interval_list_ids,

                run_clair3 = run_clair3,

                dv_threads = dv_threads,
                dv_memory = dv_memory,
                use_gpu = use_gpu,

                zones = wdl_parsable_zones
        }

        if (defined(gcs_out_dir)) {
            String smalldir = sub(select_first([gcs_out_dir]), "/$", "") + "/variants/small"
            String haptagoutdir = sub(select_first([gcs_out_dir]), "/$", "") + "/alignments"

            call FF.FinalizeToFile as FinalizeDVgVcf { input: outdir = smalldir, file = SmallVarJob.dv_g_vcf }
            call FF.FinalizeToFile as FinalizeDVgTbi { input: outdir = smalldir, file = SmallVarJob.dv_g_tbi }
            call FF.FinalizeToFile as FinalizeDVPhasedVcf { input: outdir = smalldir, file = SmallVarJob.dv_phased_vcf }
            call FF.FinalizeToFile as FinalizeDVPhasedTbi { input: outdir = smalldir, file = SmallVarJob.dv_phased_tbi }
            call FF.FinalizeToFile as FinalizeDVPhasedVcfStatusTSV { input: outdir = smalldir, file = SmallVarJob.dv_vcf_phasing_stats_tsv }
            call FF.FinalizeToFile as FinalizeDVPhasedVcfStatusGtf { input: outdir = smalldir, file = SmallVarJob.dv_vcf_phasing_stats_gtf }
            call FF.FinalizeToFile as FinalizeHapTaggedBam { input: outdir = haptagoutdir, file = SmallVarJob.haplotagged_bam }
            call FF.FinalizeToFile as FinalizeHapTaggedBai { input: outdir = haptagoutdir, file = SmallVarJob.haplotagged_bai }
            if (defined(SmallVarJob.dv_regular_resources_usage_visual)) {
                call FF.FinalizeToDir as FinalizeDVResourceUsagesVisual {
                    input: files = select_first([SmallVarJob.dv_regular_resources_usage_visual]), outdir = smalldir + "/DV_monitoring"}
            }

            if (run_clair3) {
                call FF.FinalizeToFile as FinalizeClairVcf { input: outdir = smalldir, file = select_first([SmallVarJob.clair_vcf])}
                call FF.FinalizeToFile as FinalizeClairTbi { input: outdir = smalldir, file = select_first([SmallVarJob.clair_tbi])}

                call FF.FinalizeToFile as FinalizeClairGVcf { input: outdir = smalldir, file = select_first([SmallVarJob.clair_gvcf])}
                call FF.FinalizeToFile as FinalizeClairGTbi { input: outdir = smalldir, file = select_first([SmallVarJob.clair_gtbi])}
            }
        }
    }

    ######################################################################
    # Block for SV handling
    ######################################################################
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

                zones = wdl_parsable_zones
        }

        if (defined(gcs_out_dir)) {
            String svdir = sub(select_first([gcs_out_dir]), "/$", "") + "/variants/sv"

            call FF.FinalizeToFile as FinalizePBSV { input: outdir = svdir, file = SVjob.pbsv_vcf }
            call FF.FinalizeToFile as FinalizePBSVtbi { input: outdir = svdir, file = SVjob.pbsv_tbi }

            call FF.FinalizeToFile as FinalizeSniffles { input: outdir = svdir, file = SVjob.sniffles_vcf }
            call FF.FinalizeToFile as FinalizeSnifflesTbi { input: outdir = svdir, file = SVjob.sniffles_tbi }
            call FF.FinalizeToFile as FinalizeSnifflesSnf { input: outdir = svdir, file = SVjob.sniffles_snf }
        }
    }

    ######################################################################
    # Experiment with Sniffles-2 phased SV calling
    ######################################################################
    if (call_svs && call_small_variants) {
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
        if (defined(gcs_out_dir)) {
            String svdir_copy = sub(select_first([gcs_out_dir]), "/$", "") + "/variants/sv"
            call FF.FinalizeToFile as FinalizePhasedSniffles { input: outdir = svdir_copy, file = SnifflesPhaseSV.vcf }
            call FF.FinalizeToFile as FinalizePhasedSnifflesTbi { input: outdir = svdir_copy, file = SnifflesPhaseSV.tbi }
            call FF.FinalizeToFile as FinalizePhasedSnifflesSnf { input: outdir = svdir_copy, file = SnifflesPhaseSV.snf }
        }
    }

    output {
        File? sniffles_vcf = FinalizeSniffles.gcs_path
        File? sniffles_tbi = FinalizeSnifflesTbi.gcs_path
        File? sniffles_snf = FinalizeSnifflesSnf.gcs_path

        File? sniffles_phased_vcf = FinalizePhasedSniffles.gcs_path
        File? sniffles_phased_tbi = FinalizePhasedSnifflesTbi.gcs_path
        File? sniffles_phased_snf = FinalizePhasedSnifflesSnf.gcs_path

        File? pbsv_vcf = FinalizePBSV.gcs_path
        File? pbsv_tbi = FinalizePBSVtbi.gcs_path

        File? clair_vcf = FinalizeClairVcf.gcs_path
        File? clair_tbi = FinalizeClairTbi.gcs_path
        File? clair_gvcf = FinalizeClairGVcf.gcs_path
        File? clair_gtbi = FinalizeClairGTbi.gcs_path

        File? dv_g_vcf = FinalizeDVgVcf.gcs_path
        File? dv_g_tbi = FinalizeDVgTbi.gcs_path

        File? dv_phased_vcf = FinalizeDVPhasedVcf.gcs_path
        File? dv_phased_tbi = FinalizeDVPhasedTbi.gcs_path

        File? dv_vcf_phasing_stats_tsv = FinalizeDVPhasedVcfStatusTSV.gcs_path
        File? dv_vcf_phasing_stats_gtf = FinalizeDVPhasedVcfStatusGtf.gcs_path

        File? haplotagged_bam = FinalizeHapTaggedBam.gcs_path
        File? haplotagged_bai = FinalizeHapTaggedBai.gcs_path

        String? dv_regular_resources_usage_visual = FinalizeDVResourceUsagesVisual.gcs_dir
    }
}
