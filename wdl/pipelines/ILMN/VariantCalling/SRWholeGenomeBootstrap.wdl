version 1.0

######################################################################################
## A workflow that performs single sample variant calling on Illumina reads from
## one or more flow cells. The workflow merges multiple samples into a single BAM
## prior to variant calling.
######################################################################################

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/SRUtils.wdl" as SRUTIL
import "../../../tasks/Utility/VariantUtils.wdl" as VARUTIL
import "../../../tasks/VariantCalling/CallVariantsIllumina.wdl" as VAR
import "../../../tasks/VariantCalling/HaplotypeCallerBootstrap.wdl" as HCBootstrap
import "../../../tasks/QC/AlignedMetrics.wdl" as AM
import "../../../tasks/QC/FastQC.wdl" as FastQC
import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../../tasks/QC/SampleLevelAlignedMetrics.wdl" as COV

workflow SRWholeGenome {
    input {
        Array[File] aligned_bams
        Array[File] aligned_bais

        File ref_map_file

        String participant_name

        String gcs_out_root_dir

        Boolean call_small_variants = true

        Boolean run_HC_analysis = true
        Boolean run_dv_pepper_analysis = true

        Boolean enable_hc_pileup_mode = false

        Int dvp_threads = 32
        Int dvp_memory = 128

        Int ploidy = 2

        Float heterozygosity = 0.001
        Float heterozygosity_stdev = 0.01
        Float indel_heterozygosity = 0.000125

        Int max_reads_per_alignment_start
        Int max_num_haplotypes_in_population

        File? bed_to_compute_coverage

        File? fingerprint_haplotype_db_file

        Array[String] contigs_names_to_ignore = ["RANDOM_PLACEHOLDER_VALUE"]  ## Required for ignoring any filtering - this is kind of a hack - TODO: fix the task!
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/SRWholeGenome/~{participant_name}"

    String bam_dir = outdir + "/alignments"
    String metrics_dir = outdir + "/metrics"
    String smalldir = outdir + "/variants/small"
    String recalibration_dir = outdir + "/variants/recalibration_files"

    # gather across (potential multiple) input CCS BAMs
    if (length(aligned_bams) > 1) {
        scatter (pair in zip(aligned_bams, aligned_bais)) {
            call Utils.InferSampleName as t_001_InferSampleName {input: bam = pair.left, bai = pair.right}
        }
        call Utils.CheckOnSamplenames as t_002_CheckOnSampleNames {input: sample_names = t_001_InferSampleName.sample_name}

        call Utils.MergeBams as t_003_MergeAllReads { input: bams = aligned_bams, prefix = participant_name }
    }

    File bam = select_first([t_003_MergeAllReads.merged_bam, aligned_bams[0]])
    File bai = select_first([t_003_MergeAllReads.merged_bai, aligned_bais[0]])

    # Collect sample-level metrics:
    call AM.SamStatsMap as t_004_SamStats { input: bam = bam }
    call FastQC.FastQC as t_005_FastQC { input: bam = bam, bai = bai }
    call Utils.ComputeGenomeLength as t_006_ComputeGenomeLength { input: fasta = ref_map['fasta'] }
    call SRUTIL.ComputeBamStats as t_007_ComputeBamStats { input: bam_file = bam }

    if (defined(bed_to_compute_coverage)) {
        call AM.MosDepthOverBed as t_008_MosDepth {
            input:
                bam = bam,
                bai = bai,
                bed = select_first([bed_to_compute_coverage])
        }

        call COV.SummarizeDepthOverWholeBed as t_009_RegionalCoverage {
            input:
                mosdepth_output = t_008_MosDepth.regions
        }
    }

    call FF.FinalizeToFile as t_010_FinalizeBam { input: outdir = bam_dir, file = bam, name = "~{participant_name}.bam" }
    call FF.FinalizeToFile as t_011_FinalizeBai { input: outdir = bam_dir, file = bai, name = "~{participant_name}.bam.bai" }

    if (defined(bed_to_compute_coverage)) { call FF.FinalizeToFile as t_012_FinalizeRegionalCoverage { input: outdir = bam_dir, file = select_first([t_009_RegionalCoverage.cov_summary]) } }


    call FF.FinalizeToFile as t_013_FinalizeFastQCReport {
        input:
            outdir = metrics_dir,
            file = t_005_FastQC.report
    }


    ####################################################################################################

    # Some input handling:
    if ((!run_dv_pepper_analysis) && (!run_HC_analysis)) {
        call Utils.StopWorkflow as short_variant_caller_analysis_not_provided {
            input: reason = "One of the following must be set to true: run_dv_pepper_analysis(~{run_dv_pepper_analysis}), run_HC_analysis(~{run_HC_analysis})"
        }
    }

    # Handle DeepVariant First:
    if (run_dv_pepper_analysis) {

        # Deep Variant runs better with raw base quals because it has already learned the error modes.
        # We need to revert our recalibration before calling variants:
        call SRUTIL.RevertBaseQualities as t_014_RevertBQSRQuals {
            input:
                bam = bam,
                bai = bai,
                prefix = basename(bam, ".bam") + ".reverted_base_quals"
        }

        call VAR.CallVariants as t_015_CallVariantsWithDeepVariant {
            input:
                bam               = t_014_RevertBQSRQuals.bam_out,
                bai               = t_014_RevertBQSRQuals.bai_out,
                sample_id         = participant_name,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],

                prefix = participant_name + ".deep_variant",

                call_small_variants = call_small_variants,

                run_dv_pepper_analysis = run_dv_pepper_analysis,
                dvp_threads = dvp_threads,
                dvp_memory = dvp_memory,

                mito_contig = ref_map['mt_chr_name'],
                contigs_names_to_ignore = contigs_names_to_ignore,
        }

        call FF.FinalizeToFile as t_016_FinalizeDVPepperVcf  { input: outdir = smalldir, file = select_first([t_015_CallVariantsWithDeepVariant.dvp_vcf]) }
        call FF.FinalizeToFile as t_017_FinalizeDVPepperTbi  { input: outdir = smalldir, file = select_first([t_015_CallVariantsWithDeepVariant.dvp_tbi]) }
        call FF.FinalizeToFile as t_018_FinalizeDVPepperGVcf { input: outdir = smalldir, file = select_first([t_015_CallVariantsWithDeepVariant.dvp_g_vcf]) }
        call FF.FinalizeToFile as t_019_FinalizeDVPepperGTbi { input: outdir = smalldir, file = select_first([t_015_CallVariantsWithDeepVariant.dvp_g_tbi]) }
    }

    # Now we handle HaplotypeCaller data:
    if (run_HC_analysis) {
        call HCBootstrap.CallVariantsWithHaplotypeCaller as t_020_CallVariantsWithHaplotypeCaller {
            input:
                bam               = bam,
                bai               = bai,
                sample_id         = participant_name,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                # dbsnp_vcf         = ref_map["known_sites_vcf"], #Make dbsnp_vcf optional

                ploidy            = ploidy,
                heterozygosity    = heterozygosity,
                heterozygosity_stdev = heterozygosity_stdev,
                indel_heterozygosity = indel_heterozygosity,

                max_reads_per_alignment_start = max_reads_per_alignment_start,
                max_num_haplotypes_in_population = max_num_haplotypes_in_population,

                prefix = participant_name + ".haplotype_caller",

                enable_pileup_mode = enable_hc_pileup_mode,

                mito_contig = ref_map['mt_chr_name'],
                contigs_names_to_ignore = contigs_names_to_ignore
        }

        # Make sure our sample name is correct:
        call VARUTIL.RenameSingleSampleVcf as t_021_RenameRawHcVcf {
            input:
                vcf = t_020_CallVariantsWithHaplotypeCaller.output_vcf,
                vcf_index = t_020_CallVariantsWithHaplotypeCaller.output_vcf_index,
                prefix = participant_name + ".haplotype_caller.renamed",
                new_sample_name = participant_name
        }
        call VARUTIL.RenameSingleSampleVcf as t_022_RenameRawHcGvcf {
            input:
                vcf = t_020_CallVariantsWithHaplotypeCaller.output_gvcf,
                vcf_index = t_020_CallVariantsWithHaplotypeCaller.output_gvcf_index,
                prefix = participant_name + ".haplotype_caller.renamed",
                is_gvcf = true,
                new_sample_name = participant_name
        }

        ########################################################################
        # Removed VETS
        ########################################################################
        # Fingerprinting

        if (defined(fingerprint_haplotype_db_file)) {
            call VARUTIL.ExtractFingerprintAndBarcode as t_023_FingerprintAndBarcodeVcf {
                input:
                    vcf = t_021_RenameRawHcVcf.new_sample_name_vcf,
                    vcf_index = t_021_RenameRawHcVcf.new_sample_name_vcf_index,
                    haplotype_database_file = select_first([fingerprint_haplotype_db_file]),
                    ref_fasta         = ref_map['fasta'],
                    ref_fasta_fai     = ref_map['fai'],
                    ref_dict          = ref_map['dict'],
                    prefix = participant_name
            }
        }

        ######################################################################## 
        # Removed variant filtering
        ######################################################################## 

        # Create a Keyfile for finalization:
        File keyfile = t_022_RenameRawHcGvcf.new_sample_name_vcf

        # Finalize the raw Joint Calls:
        call FF.FinalizeToFile as t_024_FinalizeHCVcf { input: outdir = smalldir, keyfile = keyfile, file = t_021_RenameRawHcVcf.new_sample_name_vcf }
        call FF.FinalizeToFile as t_025_FinalizeHCTbi { input: outdir = smalldir, keyfile = keyfile, file = t_021_RenameRawHcVcf.new_sample_name_vcf_index }
        call FF.FinalizeToFile as t_026_FinalizeHCGVcf { input: outdir = smalldir, keyfile = keyfile, file = t_022_RenameRawHcGvcf.new_sample_name_vcf }
        call FF.FinalizeToFile as t_027_FinalizeHCGTbi { input: outdir = smalldir, keyfile = keyfile, file = t_022_RenameRawHcGvcf.new_sample_name_vcf_index }
        call FF.FinalizeToFile as t_028_FinalizeHCBamOut { input: outdir = smalldir, keyfile = keyfile, file = t_020_CallVariantsWithHaplotypeCaller.bamout }
        call FF.FinalizeToFile as t_029_FinalizeHCBaiOut { input: outdir = smalldir, keyfile = keyfile, file = t_020_CallVariantsWithHaplotypeCaller.bamout_index }

        # Finalize other outputs:
        if (defined(fingerprint_haplotype_db_file)) {
            call FF.FinalizeToFile as t_030_FinalizeFingerprintVcf { input: outdir = smalldir, keyfile = keyfile, file = select_first([t_023_FingerprintAndBarcodeVcf.output_vcf]) }
        }
        
        ########################################################################
        # Removed finalization of VETS files
        ########################################################################

    }

    output {
        File aligned_bam = t_010_FinalizeBam.gcs_path
        File aligned_bai = t_011_FinalizeBai.gcs_path

        Float aligned_num_reads = t_005_FastQC.stats_map['number_of_reads']
        Float aligned_num_bases = t_004_SamStats.stats_map['bases_mapped']
        Float aligned_frac_bases = t_004_SamStats.stats_map['bases_mapped']/t_004_SamStats.stats_map['total_length']
        Float aligned_est_fold_cov = t_004_SamStats.stats_map['bases_mapped']/t_006_ComputeGenomeLength.length

        Float aligned_read_length_mean = t_005_FastQC.stats_map['read_length']

        Float insert_size_average = t_004_SamStats.stats_map['insert_size_average']
        Float insert_size_standard_deviation = t_004_SamStats.stats_map['insert_size_standard_deviation']
        Float pct_properly_paired_reads = t_004_SamStats.stats_map['percentage_of_properly_paired_reads_%']

        Float average_identity = 100.0 - (100.0*t_004_SamStats.stats_map['mismatches']/t_004_SamStats.stats_map['bases_mapped'])

        File fastqc_report = t_013_FinalizeFastQCReport.gcs_path

        Boolean successfully_processed = true

        File? bed_cov_summary = t_012_FinalizeRegionalCoverage.gcs_path

        File? fingerprint_vcf = t_030_FinalizeFingerprintVcf.gcs_path
        String? barcode = t_023_FingerprintAndBarcodeVcf.barcode
        ########################################

        File? dvp_vcf   = t_016_FinalizeDVPepperVcf.gcs_path
        File? dvp_tbi   = t_017_FinalizeDVPepperTbi.gcs_path
        File? dvp_g_vcf = t_018_FinalizeDVPepperGVcf.gcs_path
        File? dvp_g_tbi = t_019_FinalizeDVPepperGTbi.gcs_path

        ########################################
        File? hc_raw_vcf  = t_024_FinalizeHCVcf.gcs_path
        File? hc_raw_tbi  = t_025_FinalizeHCTbi.gcs_path
        File? hc_g_vcf    = t_026_FinalizeHCGVcf.gcs_path
        File? hc_g_tbi    = t_027_FinalizeHCGTbi.gcs_path
        File? hc_bamout   = t_028_FinalizeHCBamOut.gcs_path
        File? hc_baiout   = t_029_FinalizeHCBaiOut.gcs_path
    }
}
