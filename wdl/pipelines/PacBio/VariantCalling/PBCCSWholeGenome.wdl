version 1.0

import "../../../tasks/Utility/PBUtils.wdl" as PB
import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/VariantCalling/CallVariantsPBCCS.wdl" as VAR
import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/QC/SampleLevelAlignedMetrics.wdl" as COV

workflow PBCCSWholeGenome {

    meta {
        description: "A workflow that performs single sample variant calling on PacBio HiFi reads from one or more flow cells. The workflow merges multiple SMRT cells into a single BAM prior to variant calling."
    }
    parameter_meta {
        aligned_bams:       "GCS path to aligned BAM files"
        aligned_bais:       "GCS path to aligned BAM file indices"
        participant_name:   "name of the participant from whom these samples were obtained"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"

        call_svs:               "whether to call SVs"
        pbsv_discover_per_chr:  "Run the discover stage of PBSV per chromosome"

        call_small_variants: "whether to call small variants"

        run_clair3:  "to turn on Clair3 analysis or not (non-trivial increase in cost and runtime)"
        ref_scatter_interval_list_locator: "A file holding paths to interval_list files; needed only when running DV-Pepper"
        ref_scatter_interval_list_ids:     "A file that gives short IDs to the interval_list files; needed only when running DV-Pepper"

        gcp_zones: "which Google Cloud Zone to use (this has implications on how many GPUs are available and egress costs, so configure carefully)"
    }

    input {
        Array[File] aligned_bams
        Array[File] aligned_bais

        File? bed_to_compute_coverage

        File ref_map_file

        String participant_name

        String gcs_out_root_dir

        Boolean call_svs = true
        Boolean pbsv_discover_per_chr = true

        Boolean call_small_variants = true

        Boolean run_clair3 = false

        Int dv_threads = 32
        Int dv_memory = 128
        File? ref_scatter_interval_list_locator
        File? ref_scatter_interval_list_ids

        Array[String] gcp_zones = ['us-central1-a', 'us-central1-b', 'us-central1-c', 'us-central1-f']
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCCSWholeGenome/~{participant_name}"

    # gather across (potential multiple) input CCS BAMs
    if (length(aligned_bams) > 1) {
        scatter (pair in zip(aligned_bams, aligned_bais)) {
            call Utils.InferSampleName {input: bam = pair.left, bai = pair.right}
        }
        call Utils.CheckOnSamplenames {input: sample_names = InferSampleName.sample_name}

        call Utils.MergeBams as MergeAllReads { input: bams = aligned_bams, prefix = participant_name }
    }

    File bam = select_first([MergeAllReads.merged_bam, aligned_bams[0]])
    File bai = select_first([MergeAllReads.merged_bai, aligned_bais[0]])

    call COV.SampleLevelAlignedMetrics as AlignmentMetrics {
        input:
            aligned_bam = bam,
            aligned_bai = bai,
            ref_fasta   = ref_map['fasta'],
            bed_to_compute_coverage = bed_to_compute_coverage
    }

    String dir = outdir + "/alignments"

    call FF.FinalizeToFile as FinalizeBam { input: outdir = dir, file = bam, name = "~{participant_name}.bam" }
    call FF.FinalizeToFile as FinalizeBai { input: outdir = dir, file = bai, name = "~{participant_name}.bam.bai" }

    if (defined(bed_to_compute_coverage)) { call FF.FinalizeToFile as FinalizeRegionalCoverage { input: outdir = dir, file = select_first([AlignmentMetrics.bed_cov_summary]) } }

    ###########################################################
    # pacbio specific
    call PB.PBIndex as PBIndexSampleReads { input: bam = bam }
    File pbi = PBIndexSampleReads.pbi
    call FF.FinalizeToFile as FinalizePbi { input: outdir = dir, file = pbi, name = "~{participant_name}.bam.pbi" }
    ###########################################################

    ####################################################################################################
    if (call_svs || call_small_variants) {

        call GU.CollapseArrayOfStrings as get_zones {input: input_array = gcp_zones, joiner = " "}
        String wdl_parsable_zones = get_zones.collapsed

        call VAR.CallVariants {
            input:
                bam               = bam,
                bai               = bai,
                sample_id         = participant_name,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                tandem_repeat_bed = ref_map['tandem_repeat_bed'],

                prefix = participant_name,

                call_svs = call_svs,
                pbsv_discover_per_chr = pbsv_discover_per_chr,

                call_small_variants = call_small_variants,

                run_clair3 = run_clair3,

                dv_threads = dv_threads,
                dv_memory = dv_memory,

                ref_scatter_interval_list_locator = ref_scatter_interval_list_locator,
                ref_scatter_interval_list_ids = ref_scatter_interval_list_ids,

                zones = wdl_parsable_zones
        }

        String svdir = outdir + "/variants/sv"
        String smalldir = outdir + "/variants/small"

        if (call_svs) {
            call FF.FinalizeToFile as FinalizePBSV { input: outdir = svdir, file = select_first([CallVariants.pbsv_vcf]) }
            call FF.FinalizeToFile as FinalizePBSVtbi { input: outdir = svdir, file = select_first([CallVariants.pbsv_tbi]) }

            call FF.FinalizeToFile as FinalizeSniffles { input: outdir = svdir, file = select_first([CallVariants.sniffles_vcf]) }
            call FF.FinalizeToFile as FinalizeSnifflesTbi { input: outdir = svdir, file = select_first([CallVariants.sniffles_tbi]) }
        }

        if (call_small_variants) {
            if (run_clair3) {
                call FF.FinalizeToFile as FinalizeClairVcf { input: outdir = smalldir, file = select_first([CallVariants.clair_vcf])}
                call FF.FinalizeToFile as FinalizeClairTbi { input: outdir = smalldir, file = select_first([CallVariants.clair_tbi])}

                call FF.FinalizeToFile as FinalizeClairGVcf { input: outdir = smalldir, file = select_first([CallVariants.clair_gvcf])}
                call FF.FinalizeToFile as FinalizeClairGTbi { input: outdir = smalldir, file = select_first([CallVariants.clair_gtbi])}
            }

            call FF.FinalizeToFile as FinalizeDVgVcf { input: outdir = smalldir, file = select_first([CallVariants.dv_g_vcf]), name = "~{participant_name}.deepvariant.g.vcf.gz" }
            call FF.FinalizeToFile as FinalizeDVgTbi { input: outdir = smalldir, file = select_first([CallVariants.dv_g_tbi]), name = "~{participant_name}.deepvariant.g.vcf.gz.tbi" }
            call FF.FinalizeToFile as FinalizeDVPhasedVcf { input: outdir = smalldir, file = select_first([CallVariants.dv_phased_vcf]), name = "~{participant_name}.deepvariant.phased.vcf.gz" }
            call FF.FinalizeToFile as FinalizeDVPhasedTbi { input: outdir = smalldir, file = select_first([CallVariants.dv_phased_tbi]), name = "~{participant_name}.deepvariant.phased.vcf.gz.tbi" }
            call FF.FinalizeToFile as FinalizeDVPhasedVcfStatusTSV { input: outdir = smalldir, file = select_first([CallVariants.dv_vcf_phasing_stats_tsv]) }
            call FF.FinalizeToFile as FinalizeDVPhasedVcfStatusGtf { input: outdir = smalldir, file = select_first([CallVariants.dv_vcf_phasing_stats_gtf]) }
        }
    }

    call GU.GetTodayDate as today {}

    output {
        String last_preprocessing_date = today.yyyy_mm_dd

        File aligned_bam = FinalizeBam.gcs_path
        File aligned_bai = FinalizeBai.gcs_path
        File aligned_pbi = FinalizePbi.gcs_path
        Float coverage = AlignmentMetrics.coverage
        File? bed_cov_summary = FinalizeRegionalCoverage.gcs_path

        Map[String, Float] alignment_metrics = {
            'aligned_num_reads' : AlignmentMetrics.aligned_num_reads,
            'aligned_num_bases' : AlignmentMetrics.aligned_num_bases,
            'aligned_frac_bases' : AlignmentMetrics.aligned_frac_bases,
            'aligned_est_fold_cov' : AlignmentMetrics.aligned_est_fold_cov,

            'aligned_read_length_mean' : AlignmentMetrics.aligned_read_length_mean,
            'aligned_read_length_median' : AlignmentMetrics.aligned_read_length_median,
            'aligned_read_length_stdev' : AlignmentMetrics.aligned_read_length_stdev,
            'aligned_read_length_N50' : AlignmentMetrics.aligned_read_length_N50,

            'average_identity' : AlignmentMetrics.average_identity,
            'median_identity' : AlignmentMetrics.median_identity
        }

        ########################################
        File? pbsv_vcf = FinalizePBSV.gcs_path
        File? pbsv_tbi = FinalizePBSVtbi.gcs_path

        File? sniffles_vcf = FinalizeSniffles.gcs_path
        File? sniffles_tbi = FinalizeSnifflesTbi.gcs_path

        File? clair_vcf = FinalizeClairVcf.gcs_path
        File? clair_tbi = FinalizeClairTbi.gcs_path

        File? clair_gvcf = FinalizeClairGVcf.gcs_path
        File? clair_gtbi = FinalizeClairGTbi.gcs_path

        File? dv_g_vcf = FinalizeDVgVcf.gcs_path
        File? dv_g_tbi = FinalizeDVgTbi.gcs_path
        File? dv_phased_vcf = FinalizeDVPhasedVcf.gcs_path
        File? dv_phased_tbi = FinalizeDVPhasedTbi.gcs_path
        File? dv_phased_vcf_stats_tsv = FinalizeDVPhasedVcfStatusTSV.gcs_path
        File? dv_phased_vcf_stats_gtf = FinalizeDVPhasedVcfStatusGtf.gcs_path
    }
}
