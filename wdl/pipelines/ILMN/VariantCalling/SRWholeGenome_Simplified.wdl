version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/SRUtils.wdl" as SRUTIL                                          
import "../../../tasks/Utility/VariantUtils.wdl" as VARUTIL
import "../../../tasks/QC/FastQC.wdl" as FASTQC
import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../../tasks/QC/AlignedMetrics.wdl" as AM
import "../../../tasks/VariantCalling/HaplotypeCaller_Simplified.wdl" as HC
import "../../../tasks/QC/SampleLevelAlignedMetrics.wdl" as COV

workflow SRWholeGenome_Simplified {

    meta {
        author: "Jonn Smith"
        description: "This workflow performs single-sample variant calling on Illumina reads from one or more flow cells containing replicates of the same sample. The workflow merges multiple samples into a single BAM prior to variant calling."
    }

    parameter_meta {
        aligned_bams:   "Array of aligned bam files to process."
        aligned_bais:   "Array of aligned bam indices to process.  Order must correspond to `aligned_bams`."
        ref_map_file:    "Table indicating reference sequence, auxillary file locations, and metadata."
        participant_name:    "The unique identifier of this sample being processed."

        enable_hc_pileup_mode:  "If true, will enable `pileup mode` in HaplotypeCaller."
        ploidy: "Ploidy of the species being variant called."
        heterozygosity: "HaplotypeCaller Parameter - Heterozygosity value used to compute prior likelihoods for any locus. See the GATKDocs for full details on the meaning of this population genetics concept"
        heterozygosity_stdev: "HaplotypeCaller Parameter - Standard deviation of heterozygosity for SNP and indel calling."
        indel_heterozygosity: "HaplotypeCaller Parameter - Heterozygosity for indel calling. See the GATKDocs for heterozygosity for full details on the meaning of this population genetics concept"

        snp_calibration_sensitivity:    "VETS (ScoreVariantAnnotations) parameter - score below which SNP variants will be filtered."
        snp_max_unlabeled_variants: "VETS (ExtractVariantAnnotations) parameter - maximum number of unlabeled SNP variants/alleles to randomly sample with reservoir sampling.  If nonzero, annotations will also be extracted from unlabeled sites."
        snp_recalibration_annotation_values:    "VETS (ScoreSnpVariantAnnotations/ScoreVariantAnnotations) parameter - Array of annotation names to use to create the SNP variant scoring model and over which to score SNP variants."

        snp_known_reference_variants: "Array of VCF files to use as input reference variants for SNPs.  Each can be designated as either calibration or training using `snp_is_training` and `snp_is_calibration`."
        snp_known_reference_variants_index: "Array of VCF index files for `snp_known_reference_variants`.  Order should correspond to that in `snp_known_reference_variants`."
        snp_known_reference_variants_identifier: "Array of names to give to the VCF files given in `snp_known_reference_variants`.  Order should correspond to that in `snp_known_reference_variants`."
        snp_is_training: "Array of booleans indicating which files in `snp_known_reference_variants` should be used as training sets.  True -> training set.  False -> NOT a training set."
        snp_is_calibration: "Array of booleans indicating which files in `snp_known_reference_variants` should be used as calibration sets.  True ->calibration set.  False -> NOT a calibration set."

        indel_calibration_sensitivity:    "VETS (ScoreVariantAnnotations) parameter - score below which INDEL variants will be filtered."
        indel_max_unlabeled_variants: "VETS (ExtractVariantAnnotations) parameter - maximum number of unlabeled INDEL variants/alleles to randomly sample with reservoir sampling.  If nonzero, annotations will also be extracted from unlabeled sites."
        indel_recalibration_annotation_values:    "VETS (ScoreSnpVariantAnnotations/ScoreVariantAnnotations) parameter - Array of annotation names to use to create the INDEL variant scoring model and over which to score INDEL variants."

        indel_known_reference_variants: "Array of VCF files to use as input reference variants for INDELs.  Each can be designated as either calibration or training using `indel_is_training` and `indel_is_calibration`."
        indel_known_reference_variants_index: "Array of VCF index files for `indel_known_reference_variants`.  Order should correspond to that in `indel_known_reference_variants`."
        indel_known_reference_variants_identifier: "Array of names to give to the VCF files given in `indel_known_reference_variants`.  Order should correspond to that in `indel_known_reference_variants`."
        indel_is_training: "Array of booleans indicating which files in `indel_known_reference_variants` should be used as training sets.  True -> training set.  False -> NOT a training set."
        indel_is_calibration: "Array of booleans indicating which files in `indel_known_reference_variants` should be used as calibration sets.  True ->calibration set.  False -> NOT a calibration set."

        bed_to_compute_coverage: "Bed file to use as regions over which to measure coverage."

        fingerprint_haploytpe_db_file: "Haplotype DB file from which to fingerprint the input data."

        contigs_names_to_ignore:  "Array of names of contigs to ignore for the purposes of reporting variants."
        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs.  If no bucket is given, outputs will not be finalized and instead will remain in their native execution location."
    }

    input {
        Array[File] aligned_bams
        Array[File] aligned_bais

        File ref_map_file

        String participant_name

        String? gcs_out_root_dir

        Boolean enable_hc_pileup_mode = true

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

        File? interval_list

        Array[String] contigs_names_to_ignore = ["RANDOM_PLACEHOLDER_VALUE"]  ## Required for ignoring any filtering - this is kind of a hack - TODO: fix the task!
    }

    # Read ref map into map data type so we can access its fields:
    Map[String, String] ref_map = read_map(ref_map_file)

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

    # Only collect metrics if we have multiple input bam files
    if (length(aligned_bams) > 1) {
        # Collect sample-level metrics:
        call AM.SamStatsMap as SamStats { input: bam = bam }
        call FASTQC.FastQC as FastQC { input: bam = bam, bai = bai }
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

        # Calculate some final metrics that we need temporary variables for:
        Float tmp_average_identity = 100.0 - (100.0*SamStats.stats_map['mismatches']/SamStats.stats_map['bases_mapped'])
        Float tmp_aligned_frac_bases = SamStats.stats_map['bases_mapped']/SamStats.stats_map['total_length']
        Float tmp_aligned_est_fold_cov = SamStats.stats_map['bases_mapped']/ComputeGenomeLength.length
        Float tmp_aligned_num_reads = FastQC.stats_map['number_of_reads']
        Float tmp_aligned_num_bases = SamStats.stats_map['bases_mapped']
        Float tmp_aligned_read_length_mean = FastQC.stats_map['read_length']
        Float tmp_insert_size_average = SamStats.stats_map['insert_size_average']
        Float tmp_insert_size_standard_deviation = SamStats.stats_map['insert_size_standard_deviation']
        Float tmp_pct_properly_paired_reads = SamStats.stats_map['percentage_of_properly_paired_reads_%']
    }

    ####################################################################################################

    # Now we handle HaplotypeCaller data:
    call HC.CallVariantsWithHaplotypeCaller_Simplified as CallVariantsWithHaplotypeCaller {
        input:
            bam               = bam,
            bai               = bai,
            sample_id         = participant_name,
            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
            ref_dict          = ref_map['dict'],
            dbsnp_vcf         = ref_map["known_sites_vcf"],

            ploidy            = ploidy,
            heterozygosity    = heterozygosity,
            heterozygosity_stdev = heterozygosity_stdev,
            indel_heterozygosity = indel_heterozygosity,

            prefix = participant_name + ".haplotype_caller",

            enable_pileup_mode = enable_hc_pileup_mode,

            mito_contig = ref_map['mt_chr_name'],
            contigs_names_to_ignore = contigs_names_to_ignore,

            interval_list = interval_list,
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
    # Call VETS:
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

    if (defined(gcs_out_root_dir)) {

        String concrete_gcs_out_root_dir = select_first([gcs_out_root_dir])

        String outdir = sub(concrete_gcs_out_root_dir, "/$", "") + "/SRWholeGenome/~{participant_name}"

        String bam_dir = outdir + "/alignments"
        String metrics_dir = outdir + "/metrics"
        String smalldir = outdir + "/variants/small"
        String recalibration_dir = outdir + "/variants/recalibration_files"

        # Create a Keyfile for finalization:
        File keyfile = select_first([FingerprintAndBarcodeVcf.barcode_file, ScoreIndelVariantAnnotations.scored_vcf_index])

        if (length(aligned_bams) > 1) {
            call FF.FinalizeToFile as FinalizeBam { input: outdir = bam_dir, file = bam, name = "~{participant_name}.bam", keyfile=keyfile }
            call FF.FinalizeToFile as FinalizeBai { input: outdir = bam_dir, file = bai, name = "~{participant_name}.bam.bai", keyfile=keyfile }

            call FF.FinalizeToFile as FinalizeFastQCReport {
                input:
                    outdir = metrics_dir,
                    file = select_first([ FastQC.report ])
            }

            if (defined(bed_to_compute_coverage)) {
                call FF.FinalizeToFile as FinalizeRegionalCoverage { input: outdir = bam_dir, file = select_first([RegionalCoverage.cov_summary]) }
            }
        }

        # Finalize the raw Joint Calls:
        call FF.FinalizeToFile as FinalizeHCGVcf { input: outdir = smalldir, keyfile = keyfile, file = RenameRawHcGvcf.new_sample_name_vcf }
        call FF.FinalizeToFile as FinalizeHCGTbi { input: outdir = smalldir, keyfile = keyfile, file = RenameRawHcGvcf.new_sample_name_vcf_index }
        call FF.FinalizeToFile as FinalizeHCBamOut { input: outdir = smalldir, keyfile = keyfile, file = CallVariantsWithHaplotypeCaller.bamout }
        call FF.FinalizeToFile as FinalizeHCBaiOut { input: outdir = smalldir, keyfile = keyfile, file = CallVariantsWithHaplotypeCaller.bamout_index }

        # Finalize the reclibrated / filtered variants:
        call FF.FinalizeToFile as FinalizeHCRescoredVcf { input: outdir = smalldir, keyfile = keyfile, file = ScoreIndelVariantAnnotations.scored_vcf }
        call FF.FinalizeToFile as FinalizeHCRescoredTbi { input: outdir = smalldir, keyfile = keyfile, file = ScoreIndelVariantAnnotations.scored_vcf_index }

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

    # Set up output files here.
    # We have to do this because of the options for HaplotypeCaller / etc.

    # Only save the input bam files if we had to merge them:
    if (length(aligned_bams) > 1) {
        File final_aligned_bam = if defined(gcs_out_root_dir) then select_first([FinalizeBam.gcs_path]) else bam
        File final_aligned_bai = if defined(gcs_out_root_dir) then select_first([FinalizeBai.gcs_path]) else bai
        File final_fastqc_report = if defined(gcs_out_root_dir) then select_first([FinalizeFastQCReport.gcs_path]) else select_first([ FastQC.report ])

        if (defined(bed_to_compute_coverage)) {
            File final_regional_coverage_summary = if defined(gcs_out_root_dir) then select_first([FinalizeRegionalCoverage.gcs_path]) else select_first([RegionalCoverage.cov_summary])
        }
    }

    if (defined(fingerprint_haploytpe_db_file)) {
        File final_fingerprint_vcf = if defined(gcs_out_root_dir) then select_first([FinalizeFingerprintVcf.gcs_path]) else select_first([FingerprintAndBarcodeVcf.output_vcf])
    }

    output {
        ########################################
        # Primary outputs:
        File vcf = select_first([FinalizeHCRescoredVcf.gcs_path, ScoreIndelVariantAnnotations.scored_vcf])
        File tbi = select_first([FinalizeHCRescoredTbi.gcs_path, ScoreIndelVariantAnnotations.scored_vcf_index])
        File g_vcf    = select_first([FinalizeHCGVcf.gcs_path, RenameRawHcGvcf.new_sample_name_vcf])
        File g_tbi    = select_first([FinalizeHCGTbi.gcs_path, RenameRawHcGvcf.new_sample_name_vcf_index])
        File bamout   = select_first([FinalizeHCBamOut.gcs_path, CallVariantsWithHaplotypeCaller.bamout])
        File baiout   = select_first([FinalizeHCBaiOut.gcs_path, CallVariantsWithHaplotypeCaller.bamout_index])

        File? fingerprint_vcf = final_fingerprint_vcf
        String? barcode = FingerprintAndBarcodeVcf.barcode

        Boolean successfully_processed = true

        ########################################
        # Output files and statistics for multi-bam inputs:
        File? aligned_bam = final_aligned_bam
        File? aligned_bai = final_aligned_bai
        Float? average_identity = tmp_average_identity
        Float? aligned_num_reads = tmp_aligned_num_reads
        Float? aligned_num_bases = tmp_aligned_num_bases
        Float? aligned_frac_bases = tmp_aligned_frac_bases
        Float? aligned_est_fold_cov = tmp_aligned_est_fold_cov
        Float? aligned_read_length_mean = tmp_aligned_read_length_mean
        Float? insert_size_average = tmp_insert_size_average
        Float? insert_size_standard_deviation = tmp_insert_size_standard_deviation
        Float? pct_properly_paired_reads = tmp_pct_properly_paired_reads

        File? fastqc_report = final_fastqc_report
        File? bed_cov_summary = final_regional_coverage_summary
    }
}
