version 1.0

import "../../tasks/Utility/Utils.wdl" as Utils
import "../../tasks/Utility/SRUtils.wdl" as SRUTIL
import "../../tasks/Utility/Finalize.wdl" as FF
import "../../tasks/Utility/VariantUtils.wdl" as VARUTIL
import "../../tasks/Z_One_Off_Analyses/Pf_Niare_HaplotypeCaller.wdl" as Niare_HC

workflow SRWholeGenome_Pf_Niare_VETS {

    meta {
        author: "Jonn Smith"
        description: "This workflow implements a modified version of the single-sample pipeline from Niare et al. (https://doi.org/10.1186/s12936-023-04632-0) using LRMA conventions.  The modification is that this pipeline uses VETS instead of VQSR."
    }
    parameter_meta {
        aligned_bams:   "Array of aligned bam files to process."
        aligned_bais:   "Array of aligned bam indices to process.  Order must correspond to `aligned_bams`."
        ref_map_file:  "Reference map file indicating reference sequence and auxillary file locations"
        participant_name:    "The unique identifier of this sample being processed."
        vcf_calling_interval_list: "Intervals over which to call variants."
        genotype_gvcfs_intervals:  "Intervals over which to batch Joint Genotyping."

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

        contigs_names_to_ignore:  "Array of names of contigs to ignore for the purposes of reporting variants."
        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs."
    }

    input {
        Array[File] aligned_bams
        Array[File] aligned_bais

        File ref_map_file

        String participant_name

        String gcs_out_root_dir

        File vcf_calling_interval_list
        File genotype_gvcfs_intervals

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

        Array[String] contigs_names_to_ignore = ["RANDOM_PLACEHOLDER_VALUE"]  ## Required for ignoring any filtering - this is kind of a hack - TODO: fix the task!
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/SRWholeGenome_Pf_Niare_VETS/~{participant_name}"

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

    ####################################################################################################
    # HC Call Variants:

    # Now we handle HaplotypeCaller data:
    call Niare_HC.CallVariantsWithHaplotypeCaller {
        input:
            bam               = bam,
            bai               = bai,
            sample_id         = participant_name,
            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
            ref_dict          = ref_map['dict'],

            genotype_gvcfs_intervals = genotype_gvcfs_intervals,

            prefix = participant_name + ".haplotype_caller",

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

    call VARUTIL.SelectVariants as RemoveFilteredVariants {
        input:
            vcf = ScoreIndelVariantAnnotations.scored_vcf,
            vcf_index = ScoreIndelVariantAnnotations.scored_vcf_index,
            prefix = participant_name + ".filtered"
    }

    # Create a Keyfile for finalization:
    File keyfile = RemoveFilteredVariants.vcf_out_index

    # Finalize the raw Joint Calls:
    call FF.FinalizeToFile as FinalizeRawHCVcf { input: outdir = smalldir, keyfile = keyfile, file = RenameRawHcVcf.new_sample_name_vcf }
    call FF.FinalizeToFile as FinalizeRawHCTbi { input: outdir = smalldir, keyfile = keyfile, file = RenameRawHcVcf.new_sample_name_vcf_index }
    call FF.FinalizeToFile as FinalizeHCGVcf { input: outdir = smalldir, keyfile = keyfile, file = RenameRawHcGvcf.new_sample_name_vcf }
    call FF.FinalizeToFile as FinalizeHCGTbi { input: outdir = smalldir, keyfile = keyfile, file = RenameRawHcGvcf.new_sample_name_vcf_index }
    call FF.FinalizeToFile as FinalizeHCBamOut { input: outdir = smalldir, keyfile = keyfile, file = CallVariantsWithHaplotypeCaller.bamout }
    call FF.FinalizeToFile as FinalizeHCBaiOut { input: outdir = smalldir, keyfile = keyfile, file = CallVariantsWithHaplotypeCaller.bamout_index }

    # Finalize the reclibrated / filtered variants:
    call FF.FinalizeToFile as FinalizeHCRescoredVcf { input: outdir = smalldir, keyfile = keyfile, file = ScoreIndelVariantAnnotations.scored_vcf }
    call FF.FinalizeToFile as FinalizeHCRescoredTbi { input: outdir = smalldir, keyfile = keyfile, file = ScoreIndelVariantAnnotations.scored_vcf_index }

    output {
        Boolean successfully_processed = true

        ########################################

        File? hc_g_vcf    = FinalizeHCGVcf.gcs_path
        File? hc_g_tbi    = FinalizeHCGTbi.gcs_path
        File? hc_bamout   = FinalizeHCBamOut.gcs_path
        File? hc_baiout   = FinalizeHCBaiOut.gcs_path
        File? hc_raw_vcf  = FinalizeRawHCVcf.gcs_path
        File? hc_raw_tbi  = FinalizeRawHCTbi.gcs_path
        File? hc_rescored_vcf = FinalizeHCRescoredVcf.gcs_path
        File? hc_rescored_tbi = FinalizeHCRescoredVcf.gcs_path
    }
}
