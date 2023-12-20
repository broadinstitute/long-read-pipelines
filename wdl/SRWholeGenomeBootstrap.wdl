version 1.0

######################################################################################
## A workflow that performs single sample variant calling on Illumina reads from
## one or more flow cells. The workflow merges multiple samples into a single BAM
## prior to variant calling.
######################################################################################

import "tasks/Utils.wdl" as Utils
import "tasks/SRUtils.wdl" as SRUTIL
import "tasks/VariantUtils.wdl" as VARUTIL
import "tasks/CallVariantsIllumina.wdl" as VAR
import "tasks/HaplotypeCallerBootstrap.wdl" as HCBootstrap
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/FastQC.wdl" as FastQC
import "tasks/Finalize.wdl" as FF
import "tasks/SampleLevelAlignedMetrics.wdl" as COV

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

        Boolean enable_hc_pileup_mode = true

        Int dvp_threads = 32
        Int dvp_memory = 128

        Int ploidy = 2

        Float heterozygosity = 0.001
        Float heterozygosity_stdev = 0.01
        Float indel_heterozygosity = 0.000125


        Float snp_calibration_sensitivity = 0.99
        Int snp_max_unlabeled_variants = 0
        Array[String] snp_recalibration_annotation_values = [ "BaseQRankSum", "ExcessHet", "FS", "HAPCOMP", "HAPDOM", "HEC", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR", "DP" ]

        Array[File] snp_known_reference_variants
        Array[File] snp_known_reference_variants_index
        Array[File] snp_known_reference_variants_identifier
        Array[Boolean] snp_is_training
        Array[Boolean] snp_is_calibration

        Float indel_calibration_sensitivity = 0.99
        Int indel_max_unlabeled_variants = 0
        Array[String] indel_recalibration_annotation_values = [ "BaseQRankSum", "ExcessHet", "FS", "HAPCOMP", "HAPDOM", "HEC", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR", "DP" ]

        Array[File] indel_known_reference_variants
        Array[File] indel_known_reference_variants_index
        Array[File] indel_known_reference_variants_identifier
        Array[Boolean] indel_is_training
        Array[Boolean] indel_is_calibration

        File? bed_to_compute_coverage

        File? fingerprint_haploytpe_db_file

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
            call Utils.InferSampleName {input: bam = pair.left, bai = pair.right}
        }
        call Utils.CheckOnSamplenames {input: sample_names = InferSampleName.sample_name}

        call Utils.MergeBams as MergeAllReads { input: bams = aligned_bams, prefix = participant_name }
    }

    File bam = select_first([MergeAllReads.merged_bam, aligned_bams[0]])
    File bai = select_first([MergeAllReads.merged_bai, aligned_bais[0]])

    # Collect sample-level metrics:
    call AM.SamStatsMap as SamStats { input: bam = bam }
    call FastQC.FastQC as FastQC { input: bam = bam, bai = bai }
    call Utils.ComputeGenomeLength as ComputeGenomeLength { input: fasta = ref_map['fasta'] }
    call SRUTIL.ComputeBamStats as ComputeBamStats { input: bam_file = bam }

    if (defined(bed_to_compute_coverage)) {
        call AM.MosDepthOverBed as MosDepth {
            input:
                bam = bam,
                bai = bai,
                bed = select_first([bed_to_compute_coverage])
        }

        call COV.SummarizeDepthOverWholeBed as RegionalCoverage {
            input:
                mosdepth_output = MosDepth.regions
        }
    }

    call FF.FinalizeToFile as FinalizeBam { input: outdir = bam_dir, file = bam, name = "~{participant_name}.bam" }
    call FF.FinalizeToFile as FinalizeBai { input: outdir = bam_dir, file = bai, name = "~{participant_name}.bam.bai" }

    if (defined(bed_to_compute_coverage)) { call FF.FinalizeToFile as FinalizeRegionalCoverage { input: outdir = bam_dir, file = select_first([RegionalCoverage.cov_summary]) } }


    call FF.FinalizeToFile as FinalizeFastQCReport {
        input:
            outdir = metrics_dir,
            file = FastQC.report
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
        call SRUTIL.RevertBaseQualities as RevertBQSRQuals {
            input:
                bam = bam,
                bai = bai,
                prefix = basename(bam, ".bam") + ".reverted_base_quals"
        }

        call VAR.CallVariants as CallVariantsWithDeepVariant {
            input:
                bam               = RevertBQSRQuals.bam_out,
                bai               = RevertBQSRQuals.bai_out,
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

        call FF.FinalizeToFile as FinalizeDVPepperVcf  { input: outdir = smalldir, file = select_first([CallVariantsWithDeepVariant.dvp_vcf]) }
        call FF.FinalizeToFile as FinalizeDVPepperTbi  { input: outdir = smalldir, file = select_first([CallVariantsWithDeepVariant.dvp_tbi]) }
        call FF.FinalizeToFile as FinalizeDVPepperGVcf { input: outdir = smalldir, file = select_first([CallVariantsWithDeepVariant.dvp_g_vcf]) }
        call FF.FinalizeToFile as FinalizeDVPepperGTbi { input: outdir = smalldir, file = select_first([CallVariantsWithDeepVariant.dvp_g_tbi]) }
    }

    # Now we handle HaplotypeCaller data:
    if (run_HC_analysis) {
        call HCBootstrap.CallVariantsWithHaplotypeCaller {
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

                prefix = participant_name + ".haplotype_caller",

                enable_pileup_mode = enable_hc_pileup_mode,

                mito_contig = ref_map['mt_chr_name'],
                contigs_names_to_ignore = contigs_names_to_ignore,
        }

        # Make sure our sample name is correct:
        call VARUTIL.RenameSingleSampleVcf as RenameRawHcVcf {
            input:
                vcf = CallVariantsWithHaplotypeCaller.output_vcf,
                vcf_index = CallVariantsWithHaplotypeCaller.output_vcf_index,
                prefix = participant_name + ".haplotype_caller.renamed",
                new_sample_name = participant_name
        }
        call VARUTIL.RenameSingleSampleVcf as RenameRawHcGvcf {
            input:
                vcf = CallVariantsWithHaplotypeCaller.output_gvcf,
                vcf_index = CallVariantsWithHaplotypeCaller.output_gvcf_index,
                prefix = participant_name + ".haplotype_caller.renamed",
                is_gvcf = true,
                new_sample_name = participant_name
        }

        ########################################################################
        
        # Call VETS / VQSR-lite:
        call VARUTIL.ExtractVariantAnnotations as ExtractIndelVariantAnnotations {
            input:
                vcf = RenameRawHcVcf.new_sample_name_vcf,
                vcf_index = RenameRawHcVcf.new_sample_name_vcf_index,

                prefix = participant_name,
                mode = "INDEL",

                recalibration_annotation_values = indel_recalibration_annotation_values,

                known_reference_variants = indel_known_reference_variants,
                known_reference_variants_index = indel_known_reference_variants_index,
                known_reference_variants_identifier = indel_known_reference_variants_identifier,
                is_training = indel_is_training,
                is_calibration = indel_is_calibration,

                max_unlabeled_variants = indel_max_unlabeled_variants,
        }

        call VARUTIL.ExtractVariantAnnotations as ExtractSnpVariantAnnotations  {
            input:
                vcf = RenameRawHcVcf.new_sample_name_vcf,
                vcf_index = RenameRawHcVcf.new_sample_name_vcf_index,

                prefix = participant_name,
                mode = "SNP",

                recalibration_annotation_values = snp_recalibration_annotation_values,

                known_reference_variants = snp_known_reference_variants,
                known_reference_variants_index = snp_known_reference_variants_index,
                known_reference_variants_identifier = snp_known_reference_variants_identifier,
                is_training = snp_is_training,
                is_calibration = snp_is_calibration,

                max_unlabeled_variants = snp_max_unlabeled_variants,
        }

        call VARUTIL.TrainVariantAnnotationsModel as TrainIndelVariantAnnotationsModel {
            input:
                annotation_hdf5 = ExtractIndelVariantAnnotations.annotation_hdf5,
                mode = "INDEL",
                prefix = participant_name,
        }

        call VARUTIL.TrainVariantAnnotationsModel as TrainSnpVariantAnnotationsModel {
            input:
                annotation_hdf5 = ExtractSnpVariantAnnotations.annotation_hdf5,
                mode = "SNP",
                prefix = participant_name,
        }

        call VARUTIL.ScoreVariantAnnotations as ScoreSnpVariantAnnotations {
            input:
                vcf = RenameRawHcVcf.new_sample_name_vcf,
                vcf_index = RenameRawHcVcf.new_sample_name_vcf_index,

                sites_only_extracted_vcf = ExtractSnpVariantAnnotations.sites_only_vcf,
                sites_only_extracted_vcf_index = ExtractSnpVariantAnnotations.sites_only_vcf_index,

                model_prefix = participant_name + "_train_SNP",
                model_files = flatten([[TrainSnpVariantAnnotationsModel.training_scores, TrainSnpVariantAnnotationsModel.positive_model_scorer_pickle], select_all([
                    TrainSnpVariantAnnotationsModel.unlabeled_positive_model_scores,
                    TrainSnpVariantAnnotationsModel.calibration_set_scores,
                    TrainSnpVariantAnnotationsModel.negative_model_scorer_pickle
                ])]),
                prefix = participant_name + "_SNP",
                mode = "SNP",

                calibration_sensitivity_threshold = snp_calibration_sensitivity,

                recalibration_annotation_values = snp_recalibration_annotation_values,

                known_reference_variants = snp_known_reference_variants,
                known_reference_variants_index = snp_known_reference_variants_index,
                known_reference_variants_identifier = snp_known_reference_variants_identifier,
                is_training = snp_is_training,
                is_calibration = snp_is_calibration,
        }

        call VARUTIL.ScoreVariantAnnotations as ScoreIndelVariantAnnotations {
            input:
                vcf = ScoreSnpVariantAnnotations.scored_vcf,
                vcf_index = ScoreSnpVariantAnnotations.scored_vcf_index,

                sites_only_extracted_vcf = ExtractIndelVariantAnnotations.sites_only_vcf,
                sites_only_extracted_vcf_index = ExtractIndelVariantAnnotations.sites_only_vcf_index,

                model_prefix = participant_name + "_train_INDEL",
                model_files = flatten([[TrainIndelVariantAnnotationsModel.training_scores, TrainIndelVariantAnnotationsModel.positive_model_scorer_pickle], select_all([
                    TrainIndelVariantAnnotationsModel.unlabeled_positive_model_scores,
                    TrainIndelVariantAnnotationsModel.calibration_set_scores,
                    TrainIndelVariantAnnotationsModel.negative_model_scorer_pickle
                ])]),
                prefix = participant_name + "_ALL",
                mode = "INDEL",

                calibration_sensitivity_threshold = indel_calibration_sensitivity,

                recalibration_annotation_values = indel_recalibration_annotation_values,

                known_reference_variants = indel_known_reference_variants,
                known_reference_variants_index = indel_known_reference_variants_index,
                known_reference_variants_identifier = indel_known_reference_variants_identifier,
                is_training = indel_is_training,
                is_calibration = indel_is_calibration,
        }
        ########################################################################

        if (defined(fingerprint_haploytpe_db_file)) {
            call VARUTIL.ExtractFingerprintAndBarcode as FingerprintAndBarcodeVcf {
                input:
                    vcf = ScoreIndelVariantAnnotations.scored_vcf,
                    vcf_index = ScoreIndelVariantAnnotations.scored_vcf_index,
                    haplotype_database_file = select_first([fingerprint_haploytpe_db_file]),
                    ref_fasta         = ref_map['fasta'],
                    ref_fasta_fai     = ref_map['fai'],
                    ref_dict          = ref_map['dict'],
                    prefix = participant_name
            }
        }

        call VARUTIL.SelectVariants as RemoveFilteredVariants {
            input:
                vcf = ScoreIndelVariantAnnotations.scored_vcf,
                vcf_index = ScoreIndelVariantAnnotations.scored_vcf_index,
                prefix = participant_name + ".filtered"
        }

        # Create a Keyfile for finalization:
        File keyfile = RemoveFilteredVariants.vcf_out_index

        # Finalize the raw Joint Calls:
        call FF.FinalizeToFile as FinalizeHCVcf { input: outdir = smalldir, keyfile = keyfile, file = RenameRawHcVcf.new_sample_name_vcf }
        call FF.FinalizeToFile as FinalizeHCTbi { input: outdir = smalldir, keyfile = keyfile, file = RenameRawHcVcf.new_sample_name_vcf_index }
        call FF.FinalizeToFile as FinalizeHCGVcf { input: outdir = smalldir, keyfile = keyfile, file = RenameRawHcGvcf.new_sample_name_vcf }
        call FF.FinalizeToFile as FinalizeHCGTbi { input: outdir = smalldir, keyfile = keyfile, file = RenameRawHcGvcf.new_sample_name_vcf_index }
        call FF.FinalizeToFile as FinalizeHCBamOut { input: outdir = smalldir, keyfile = keyfile, file = CallVariantsWithHaplotypeCaller.bamout }
        call FF.FinalizeToFile as FinalizeHCBaiOut { input: outdir = smalldir, keyfile = keyfile, file = CallVariantsWithHaplotypeCaller.bamout_index }

        # Finalize the reclibrated / filtered variants:
        call FF.FinalizeToFile as FinalizeHCRescoredVcf { input: outdir = smalldir, keyfile = keyfile, file = ScoreIndelVariantAnnotations.scored_vcf }
        call FF.FinalizeToFile as FinalizeHCRescoredTbi { input: outdir = smalldir, keyfile = keyfile, file = ScoreIndelVariantAnnotations.scored_vcf_index }
        call FF.FinalizeToFile as FinalizeHCRescoredFilteredVcf { input: outdir = smalldir, keyfile = keyfile, file = RemoveFilteredVariants.vcf_out }
        call FF.FinalizeToFile as FinalizeHCRescoredFilteredTbi { input: outdir = smalldir, keyfile = keyfile, file = RemoveFilteredVariants.vcf_out_index }

        # Finalize other outputs:
        if (defined(fingerprint_haploytpe_db_file)) {
            call FF.FinalizeToFile as FinalizeFingerprintVcf { input: outdir = smalldir, keyfile = keyfile, file = select_first([FingerprintAndBarcodeVcf.output_vcf]) }
        }

        ################################
        # Finalize the VETS files:
        ############

        # ExtractVariantAnnotations:
        call FF.FinalizeToFile as FinalizeSnpExtractedAnnotations { input: outdir = recalibration_dir, keyfile = keyfile, file = ExtractSnpVariantAnnotations.annotation_hdf5 }
        call FF.FinalizeToFile as FinalizeSnpExtractedSitesOnlyVcf { input: outdir = recalibration_dir, keyfile = keyfile, file = ExtractSnpVariantAnnotations.sites_only_vcf }
        call FF.FinalizeToFile as FinalizeSnpExtractedSitesOnlyVcfIndex { input: outdir = recalibration_dir, keyfile = keyfile, file = ExtractSnpVariantAnnotations.sites_only_vcf_index }
        if (defined(ExtractSnpVariantAnnotations.unlabeled_annotation_hdf5)) {
            call FF.FinalizeToFile as FinalizeSnpExtractedUnlabeledAnnotations { input: outdir = recalibration_dir, keyfile = keyfile, file = select_first([ExtractSnpVariantAnnotations.unlabeled_annotation_hdf5]) }
        }
        call FF.FinalizeToFile as FinalizeIndelExtractedAnnotations { input: outdir = recalibration_dir, keyfile = keyfile, file = ExtractIndelVariantAnnotations.annotation_hdf5 }
        call FF.FinalizeToFile as FinalizeIndelExtractedSitesOnlyVcf { input: outdir = recalibration_dir, keyfile = keyfile, file = ExtractIndelVariantAnnotations.sites_only_vcf }
        call FF.FinalizeToFile as FinalizeIndelExtractedSitesOnlyVcfIndex { input: outdir = recalibration_dir, keyfile = keyfile, file = ExtractIndelVariantAnnotations.sites_only_vcf_index }
        if (defined(ExtractIndelVariantAnnotations.unlabeled_annotation_hdf5)) {
            call FF.FinalizeToFile as FinalizeIndelExtractedUnlabeledAnnotations { input: outdir = recalibration_dir, keyfile = keyfile, file = select_first([ExtractIndelVariantAnnotations.unlabeled_annotation_hdf5]) }
        }

        # TrainVariantAnnotationsModel
        call FF.FinalizeToFile as FinalizeSnpTrainVariantAnnotationsTrainingScores { input: outdir = recalibration_dir, keyfile = keyfile, file = TrainSnpVariantAnnotationsModel.training_scores }
        call FF.FinalizeToFile as FinalizeSnpTrainVariantAnnotationsPositiveModelScorer { input: outdir = recalibration_dir, keyfile = keyfile, file = TrainSnpVariantAnnotationsModel.positive_model_scorer_pickle }
        if (defined(TrainSnpVariantAnnotationsModel.unlabeled_positive_model_scores)) {
            call FF.FinalizeToFile as FinalizeSnpTrainVariantAnnotationsUnlabeledPositiveModelScores { input: outdir = recalibration_dir, keyfile = keyfile, file = select_first([TrainSnpVariantAnnotationsModel.unlabeled_positive_model_scores]) }
        }
        if (defined(TrainSnpVariantAnnotationsModel.calibration_set_scores)) {
            call FF.FinalizeToFile as FinalizeSnpTrainVariantAnnotationsCalibrationSetScores { input: outdir = recalibration_dir, keyfile = keyfile, file = select_first([TrainSnpVariantAnnotationsModel.calibration_set_scores]) }
        }
        if (defined(TrainSnpVariantAnnotationsModel.negative_model_scorer_pickle)) {
            call FF.FinalizeToFile as FinalizeSnpTrainVariantAnnotationsNegativeModelScorer { input: outdir = recalibration_dir, keyfile = keyfile, file = select_first([TrainSnpVariantAnnotationsModel.negative_model_scorer_pickle]) }
        }

        call FF.FinalizeToFile as FinalizeIndelTrainVariantAnnotationsTrainingScores { input: outdir = recalibration_dir, keyfile = keyfile, file = TrainIndelVariantAnnotationsModel.training_scores }
        call FF.FinalizeToFile as FinalizeIndelTrainVariantAnnotationsPositiveModelScorer { input: outdir = recalibration_dir, keyfile = keyfile, file = TrainIndelVariantAnnotationsModel.positive_model_scorer_pickle }
        if (defined(TrainIndelVariantAnnotationsModel.unlabeled_positive_model_scores)) {
            call FF.FinalizeToFile as FinalizeIndelTrainVariantAnnotationsUnlabeledPositiveModelScores { input: outdir = recalibration_dir, keyfile = keyfile, file = select_first([TrainIndelVariantAnnotationsModel.unlabeled_positive_model_scores]) }
        }
        if (defined(TrainIndelVariantAnnotationsModel.calibration_set_scores)) {
            call FF.FinalizeToFile as FinalizeIndelTrainVariantAnnotationsCalibrationSetScores { input: outdir = recalibration_dir, keyfile = keyfile, file = select_first([TrainIndelVariantAnnotationsModel.calibration_set_scores]) }
        }
        if (defined(TrainIndelVariantAnnotationsModel.negative_model_scorer_pickle)) {
            call FF.FinalizeToFile as FinalizeIndelTrainVariantAnnotationsNegativeModelScorer { input: outdir = recalibration_dir, keyfile = keyfile, file = select_first([TrainIndelVariantAnnotationsModel.negative_model_scorer_pickle]) }
        }

        # ScoreVariantAnnotations
        call FF.FinalizeToFile as FinalizeScoreSnpVariantAnnotationsScoredVcf { input: outdir = recalibration_dir, keyfile = keyfile, file = ScoreSnpVariantAnnotations.scored_vcf }
        call FF.FinalizeToFile as FinalizeScoreSnpVariantAnnotationsScoredVcfIndex { input: outdir = recalibration_dir, keyfile = keyfile, file = ScoreSnpVariantAnnotations.scored_vcf_index }
        if (defined(ScoreSnpVariantAnnotations.annotations_hdf5)) {
            call FF.FinalizeToFile as FinalizeScoreSnpVariantAnnotationsAnnotationsHdf5 { input: outdir = recalibration_dir, keyfile = keyfile, file = select_first([ScoreSnpVariantAnnotations.annotations_hdf5]) }
        }
        if (defined(ScoreSnpVariantAnnotations.scores_hdf5)) {
            call FF.FinalizeToFile as FinalizeScoreSnpVariantAnnotationsScoresHdf5 { input: outdir = recalibration_dir, keyfile = keyfile, file = select_first([ScoreSnpVariantAnnotations.scores_hdf5]) }
        }

        call FF.FinalizeToFile as FinalizeScoreIndelVariantAnnotationsScoredVcf { input: outdir = recalibration_dir, keyfile = keyfile, file = ScoreIndelVariantAnnotations.scored_vcf }
        call FF.FinalizeToFile as FinalizeScoreIndelVariantAnnotationsScoredVcfIndex { input: outdir = recalibration_dir, keyfile = keyfile, file = ScoreIndelVariantAnnotations.scored_vcf_index }
        if (defined(ScoreIndelVariantAnnotations.annotations_hdf5)) {
            call FF.FinalizeToFile as FinalizeScoreIndelVariantAnnotationsAnnotationsHdf5 { input: outdir = recalibration_dir, keyfile = keyfile, file = select_first([ScoreIndelVariantAnnotations.annotations_hdf5]) }
        }
        if (defined(ScoreIndelVariantAnnotations.scores_hdf5)) {
            call FF.FinalizeToFile as FinalizeScoreIndelVariantAnnotationsScoresHdf5 { input: outdir = recalibration_dir, keyfile = keyfile, file = select_first([ScoreIndelVariantAnnotations.scores_hdf5]) }
        }
    }

    output {
        File aligned_bam = FinalizeBam.gcs_path
        File aligned_bai = FinalizeBai.gcs_path

        Float aligned_num_reads = FastQC.stats_map['number_of_reads']
        Float aligned_num_bases = SamStats.stats_map['bases_mapped']
        Float aligned_frac_bases = SamStats.stats_map['bases_mapped']/SamStats.stats_map['total_length']
        Float aligned_est_fold_cov = SamStats.stats_map['bases_mapped']/ComputeGenomeLength.length

        Float aligned_read_length_mean = FastQC.stats_map['read_length']

        Float insert_size_average = SamStats.stats_map['insert_size_average']
        Float insert_size_standard_deviation = SamStats.stats_map['insert_size_standard_deviation']
        Float pct_properly_paired_reads = SamStats.stats_map['percentage_of_properly_paired_reads_%']

        Float average_identity = 100.0 - (100.0*SamStats.stats_map['mismatches']/SamStats.stats_map['bases_mapped'])

        File fastqc_report = FinalizeFastQCReport.gcs_path

        Boolean successfully_processed = true

        File? bed_cov_summary = FinalizeRegionalCoverage.gcs_path

        File? fingerprint_vcf = FinalizeFingerprintVcf.gcs_path
        String? barcode = FingerprintAndBarcodeVcf.barcode

        ########################################

        File? dvp_vcf   = FinalizeDVPepperVcf.gcs_path
        File? dvp_tbi   = FinalizeDVPepperTbi.gcs_path
        File? dvp_g_vcf = FinalizeDVPepperGVcf.gcs_path
        File? dvp_g_tbi = FinalizeDVPepperGTbi.gcs_path

        ########################################

        File? hc_g_vcf    = FinalizeHCGVcf.gcs_path
        File? hc_g_tbi    = FinalizeHCGTbi.gcs_path
        File? hc_bamout   = FinalizeHCBamOut.gcs_path
        File? hc_baiout   = FinalizeHCBaiOut.gcs_path
        File? hc_raw_vcf  = FinalizeHCVcf.gcs_path
        File? hc_raw_tbi  = FinalizeHCTbi.gcs_path
        File? hc_rescored_vcf = FinalizeHCRescoredVcf.gcs_path
        File? hc_rescored_tbi = FinalizeHCRescoredTbi.gcs_path
        File? hc_rescored_filtered_vcf = FinalizeHCRescoredFilteredVcf.gcs_path
        File? hc_rescored_filtered_tbi = FinalizeHCRescoredFilteredTbi.gcs_path
    }
}
