version 1.0

#############################################################################################################
## A workflow that performs joint calling on single-sample gVCFs from GATK4 HaplotypeCaller using GenomicsDB.
#############################################################################################################

import "tasks/SRJointGenotyping.wdl" as SRJOINT
import "tasks/VariantUtils.wdl" as VARUTIL
import "tasks/Utils.wdl" as UTILS
import "tasks/Hail.wdl" as Hail
import "tasks/FunctionalAnnotation.wdl" as FUNK
import "tasks/SGKit.wdl" as SGKit
import "tasks/Finalize.wdl" as FF

workflow SRJointCallGVCFsWithGenomicsDB {
    input {
        Array[File] gvcfs
        Array[File] gvcf_indices

        File ref_map_file

        File interval_list

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

        Array[File]?   annotation_bed_files
        Array[File]?   annotation_bed_file_indexes
        Array[String]? annotation_bed_file_annotation_names

        File? snpeff_db

        String prefix

        String gcs_out_root_dir
    }

    parameter_meta {
        gvcfs:            "GCS paths to gVCF files"
        gvcf_indices:     "GCS paths to gVCF tbi files"
        ref_map_file:     "table indicating reference sequence and auxillary file locations"
        prefix:           "prefix for output joint-called gVCF and tabix index"
        gcs_out_root_dir: "GCS bucket to store the reads, variants, and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/SRJointCallGVCFsWithGenomicsDB/~{prefix}"

    Map[String, String] ref_map = read_map(ref_map_file)

    # Create interval list over which to shard the processing:
    call UTILS.MakeChrIntervalList as MakeChrIntervalList {
        input:
            ref_dict = ref_map['dict'],
    }

    # Create sample-name map:
    call SRJOINT.CreateSampleNameMap as CreateSampleNameMap {
        input:
            gvcfs = gvcfs,
            prefix = prefix
    }

    # Shard by contig for speed:
    scatter (idx_1 in range(length(MakeChrIntervalList.contig_interval_list_files))) {

        String contig = MakeChrIntervalList.chrs[idx_1][0]
        File contig_interval_list = MakeChrIntervalList.contig_interval_list_files[idx_1]

        # Import our data into GenomicsDB:
        call SRJOINT.ImportGVCFs as ImportGVCFsIntoGenomicsDB {
            input:
                sample_name_map = CreateSampleNameMap.sample_name_map,
                interval_list   = contig_interval_list,
                ref_fasta       = ref_map['fasta'],
                ref_fasta_fai   = ref_map['fai'],
                ref_dict        = ref_map['dict'],
                prefix          = prefix + "." + contig,
                batch_size      = 50,
                # We need to override this because we're not actually sending the GVCF over (just a list)
                # ALSO, we're currently tarring the genomicsDB, so we need at least double the space here, plus some slop:
                runtime_attr_override = object {disk_gb: 10 + (3 * CreateSampleNameMap.total_gvcf_size_gb) + (2 * ceil(size(ref_map['fasta'], "GB"))), preemptible_tries: 0}
        }

        # Joint call
        call SRJOINT.GenotypeGVCFs as JointCallGVCFs {
            input:
                input_gvcf_data = ImportGVCFsIntoGenomicsDB.output_genomicsdb,
                interval_list   = contig_interval_list,
                ref_fasta       = ref_map['fasta'],
                ref_fasta_fai   = ref_map['fai'],
                ref_dict        = ref_map['dict'],
                dbsnp_vcf       = ref_map["known_sites_vcf"],
                prefix          = prefix + "." + contig + ".raw",
                runtime_attr_override = object {preemptible_tries: 0},  # Disable preemption for prototype.
        }

        # First make a sites-only VCF for recal (smaller file, easier to work with):
        call VARUTIL.MakeSitesOnlyVcf as MakeSitesOnlyVCF {
            input:
                vcf = JointCallGVCFs.output_vcf,
                vcf_index = JointCallGVCFs.output_vcf_index,
                prefix = prefix + "." + contig + ".sites_only"
        }
    }

    # Merge all sites-only VCFs
    call VARUTIL.GatherVcfs as MergeSitesOnlyVCFs {
        input:
            input_vcfs = MakeSitesOnlyVCF.sites_only_vcf,
            input_vcf_indices = MakeSitesOnlyVCF.sites_only_vcf_index,
            prefix = prefix + ".sites_only"
    }

    ########################################################################
    # Call VETS / VQSR-lite:
    call VARUTIL.ExtractVariantAnnotations as ExtractIndelVariantAnnotations {
        input:
            vcf = MergeSitesOnlyVCFs.output_vcf,
            vcf_index = MergeSitesOnlyVCFs.output_vcf_index,

            prefix = prefix,
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
            vcf = MergeSitesOnlyVCFs.output_vcf,
            vcf_index = MergeSitesOnlyVCFs.output_vcf_index,

            prefix = prefix,
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
            prefix = prefix,
    }

    call VARUTIL.TrainVariantAnnotationsModel as TrainSnpVariantAnnotationsModel {
        input:
            annotation_hdf5 = ExtractIndelVariantAnnotations.annotation_hdf5,
            mode = "INDEL",
            prefix = prefix,
    }

    # Shard by contig for speed:
    scatter (idx_2 in range(length(JointCallGVCFs.output_vcf))) {

        String contig_2 = MakeChrIntervalList.chrs[idx_2][0]
        File joint_called_vcf = JointCallGVCFs.output_vcf[idx_2]
        File joint_called_vcf_index = JointCallGVCFs.output_vcf_index[idx_2]

        call VARUTIL.ScoreVariantAnnotations as ScoreSnpVariantAnnotations {
            input:
                vcf = joint_called_vcf,
                vcf_index = joint_called_vcf_index,

                sites_only_extracted_vcf = ExtractSnpVariantAnnotations.sites_only_vcf,
                sites_only_extracted_vcf_index = ExtractSnpVariantAnnotations.sites_only_vcf_index,

                model_prefix = prefix + "_train_SNP",
                model_files = flatten([[TrainSnpVariantAnnotationsModel.training_scores, TrainSnpVariantAnnotationsModel.positive_model_scorer_pickle], select_all([
                    TrainSnpVariantAnnotationsModel.unlabeled_positive_model_scores,
                    TrainSnpVariantAnnotationsModel.calibration_set_scores,
                    TrainSnpVariantAnnotationsModel.negative_model_scorer_pickle
                ])]),
                prefix = prefix + "_SNP_" + contig_2,
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

                model_prefix = prefix + "_train_INDEL",
                model_files = flatten([[TrainIndelVariantAnnotationsModel.training_scores, TrainIndelVariantAnnotationsModel.positive_model_scorer_pickle], select_all([
                    TrainIndelVariantAnnotationsModel.unlabeled_positive_model_scores,
                    TrainIndelVariantAnnotationsModel.calibration_set_scores,
                    TrainIndelVariantAnnotationsModel.negative_model_scorer_pickle
                ])]),
                prefix = prefix + "_ALL_" + contig_2,
                mode = "INDEL",

                calibration_sensitivity_threshold = indel_calibration_sensitivity,

                recalibration_annotation_values = indel_recalibration_annotation_values,

                known_reference_variants = indel_known_reference_variants,
                known_reference_variants_index = indel_known_reference_variants_index,
                known_reference_variants_identifier = indel_known_reference_variants_identifier,
                is_training = indel_is_training,
                is_calibration = indel_is_calibration,
        }

        # Now we need to annotate our variants by region:
        if (defined(annotation_bed_files)) {
            call VARUTIL.AnnotateVcfWithBedRegions as AnnotateVcfRegions {
                input:
                    vcf = ScoreIndelVariantAnnotations.scored_vcf,
                    vcf_index = ScoreIndelVariantAnnotations.scored_vcf_index,
                    bed_files = select_first([annotation_bed_files]),
                    bed_file_indexes = select_first([annotation_bed_file_indexes]),
                    bed_file_annotation_names = select_first([annotation_bed_file_annotation_names]),
                    prefix = basename(basename(ScoreIndelVariantAnnotations.scored_vcf, ".vcf.gz"), ".vcf") + ".region_annotated",
            }
        }

        File recalibrated_vcf = select_first([AnnotateVcfRegions.annotated_vcf, ScoreIndelVariantAnnotations.scored_vcf])
        File recalibrated_vcf_index = select_first([AnnotateVcfRegions.annotated_vcf_index, ScoreIndelVariantAnnotations.scored_vcf_index])

        # Now functionally annotate each VCF:
        if (defined(snpeff_db)) {
            call FUNK.FunctionallyAnnotateVariants as FunctionallyAnnotate {
                input:
                    vcf = recalibrated_vcf,
                    snpeff_db = select_first([snpeff_db])
            }
            call VARUTIL.IndexVCF as IndexFunkyVcf {
                input:
                    vcf = FunctionallyAnnotate.annotated_vcf
            }
        }

        File vcf_for_merging = select_first([FunctionallyAnnotate.annotated_vcf, recalibrated_vcf])
        File vcf_index_for_merging = select_first([IndexFunkyVcf.tbi, recalibrated_vcf_index])
    }

    # Consolidate files:
    call VARUTIL.GatherVcfs as GatherRawVcfs {
        input:
            input_vcfs = JointCallGVCFs.output_vcf,
            input_vcf_indices = JointCallGVCFs.output_vcf_index,
            prefix = prefix + ".raw.combined"
    }

    # Consolidate files:
    call VARUTIL.GatherVcfs as GatherRescoredVcfs {
        input:
            input_vcfs = vcf_for_merging,
            input_vcf_indices = vcf_index_for_merging,
            prefix = prefix + ".rescored.combined"
    }

    # Convert to Zarr
    call SGKit.ConvertToZarrStore as ConvertToZarr {
        input:
            gvcf = GatherRescoredVcfs.output_vcf,
            tbi = GatherRescoredVcfs.output_vcf_index,
            prefix = prefix,
            outdir = outdir
    }

    # Convert the output to a HAIL Matrix Table:
    call Hail.ConvertToHailMT as CreateHailMatrixTable {
        input:
            gvcf = GatherRescoredVcfs.output_vcf,
            tbi = GatherRescoredVcfs.output_vcf_index,
            reference = sub(sub(ref_map["fasta"], "^.*/", ""), "\.[fasta]*$", ""),
            ref_fasta = ref_map["fasta"],
            ref_fai = ref_map["fai"],
            prefix = prefix,
            outdir = outdir
    }

    ################################
    # Finalize the regular output files:
    ############
    File keyfile = CreateHailMatrixTable.monitoring_log
    String recalibration_dir = outdir + "/recalibration_files"
    String recalibration_model_dir = outdir + "/recalibration_files/model"
    String recalibration_results_dir = outdir + "/recalibration_files/results"

    call FF.FinalizeToDir as FinalizeGenomicsDB { input: outdir = outdir + "/GenomicsDB", keyfile = keyfile, files = ImportGVCFsIntoGenomicsDB.output_genomicsdb }

    call FF.FinalizeToFile as FinalizeRawVCF { input: outdir = outdir, keyfile = keyfile, file = GatherRawVcfs.output_vcf }
    call FF.FinalizeToFile as FinalizeRawTBI { input: outdir = outdir, keyfile = keyfile, file = GatherRawVcfs.output_vcf_index }

    call FF.FinalizeToFile as FinalizeVETSVCF { input: outdir = outdir, keyfile = keyfile, file = GatherRescoredVcfs.output_vcf }
    call FF.FinalizeToFile as FinalizeVETSTBI { input: outdir = outdir, keyfile = keyfile, file = GatherRescoredVcfs.output_vcf_index }

    ################################
    # Finalize the VETS files:
    ############

    # ExtractVariantAnnotations:
    call FF.FinalizeToFile as FinalizeSnpExtractedAnnotations { input: outdir = recalibration_model_dir, keyfile = keyfile, file = ExtractSnpVariantAnnotations.annotation_hdf5 }
    call FF.FinalizeToFile as FinalizeSnpExtractedSitesOnlyVcf { input: outdir = recalibration_model_dir, keyfile = keyfile, file = ExtractSnpVariantAnnotations.sites_only_vcf }
    call FF.FinalizeToFile as FinalizeSnpExtractedSitesOnlyVcfIndex { input: outdir = recalibration_model_dir, keyfile = keyfile, file = ExtractSnpVariantAnnotations.sites_only_vcf_index }
    if (defined(ExtractSnpVariantAnnotations.unlabeled_annotation_hdf5)) {
        call FF.FinalizeToFile as FinalizeSnpExtractedUnlabeledAnnotations { input: outdir = recalibration_model_dir, keyfile = keyfile, file = select_first([ExtractSnpVariantAnnotations.unlabeled_annotation_hdf5]) }
    }
    call FF.FinalizeToFile as FinalizeIndelExtractedAnnotations { input: outdir = recalibration_model_dir, keyfile = keyfile, file = ExtractIndelVariantAnnotations.annotation_hdf5 }
    call FF.FinalizeToFile as FinalizeIndelExtractedSitesOnlyVcf { input: outdir = recalibration_model_dir, keyfile = keyfile, file = ExtractIndelVariantAnnotations.sites_only_vcf }
    call FF.FinalizeToFile as FinalizeIndelExtractedSitesOnlyVcfIndex { input: outdir = recalibration_model_dir, keyfile = keyfile, file = ExtractIndelVariantAnnotations.sites_only_vcf_index }
    if (defined(ExtractIndelVariantAnnotations.unlabeled_annotation_hdf5)) {
        call FF.FinalizeToFile as FinalizeIndelExtractedUnlabeledAnnotations { input: outdir = recalibration_model_dir, keyfile = keyfile, file = select_first([ExtractIndelVariantAnnotations.unlabeled_annotation_hdf5]) }
    }

    # TrainVariantAnnotationsModel
    call FF.FinalizeToFile as FinalizeSnpTrainVariantAnnotationsTrainingScores { input: outdir = recalibration_model_dir, keyfile = keyfile, file = TrainSnpVariantAnnotationsModel.training_scores }
    call FF.FinalizeToFile as FinalizeSnpTrainVariantAnnotationsPositiveModelScorer { input: outdir = recalibration_model_dir, keyfile = keyfile, file = TrainSnpVariantAnnotationsModel.positive_model_scorer_pickle }
    if (defined(TrainSnpVariantAnnotationsModel.unlabeled_positive_model_scores)) {
        call FF.FinalizeToFile as FinalizeSnpTrainVariantAnnotationsUnlabeledPositiveModelScores { input: outdir = recalibration_model_dir, keyfile = keyfile, file = select_first([TrainSnpVariantAnnotationsModel.unlabeled_positive_model_scores]) }
    }
    if (defined(TrainSnpVariantAnnotationsModel.calibration_set_scores)) {
        call FF.FinalizeToFile as FinalizeSnpTrainVariantAnnotationsCalibrationSetScores { input: outdir = recalibration_model_dir, keyfile = keyfile, file = select_first([TrainSnpVariantAnnotationsModel.calibration_set_scores]) }
    }
    if (defined(TrainSnpVariantAnnotationsModel.negative_model_scorer_pickle)) {
        call FF.FinalizeToFile as FinalizeSnpTrainVariantAnnotationsNegativeModelScorer { input: outdir = recalibration_model_dir, keyfile = keyfile, file = select_first([TrainSnpVariantAnnotationsModel.negative_model_scorer_pickle]) }
    }

    call FF.FinalizeToFile as FinalizeIndelTrainVariantAnnotationsTrainingScores { input: outdir = recalibration_model_dir, keyfile = keyfile, file = TrainIndelVariantAnnotationsModel.training_scores }
    call FF.FinalizeToFile as FinalizeIndelTrainVariantAnnotationsPositiveModelScorer { input: outdir = recalibration_model_dir, keyfile = keyfile, file = TrainIndelVariantAnnotationsModel.positive_model_scorer_pickle }
    if (defined(TrainIndelVariantAnnotationsModel.unlabeled_positive_model_scores)) {
        call FF.FinalizeToFile as FinalizeIndelTrainVariantAnnotationsUnlabeledPositiveModelScores { input: outdir = recalibration_model_dir, keyfile = keyfile, file = select_first([TrainIndelVariantAnnotationsModel.unlabeled_positive_model_scores]) }
    }
    if (defined(TrainIndelVariantAnnotationsModel.calibration_set_scores)) {
        call FF.FinalizeToFile as FinalizeIndelTrainVariantAnnotationsCalibrationSetScores { input: outdir = recalibration_model_dir, keyfile = keyfile, file = select_first([TrainIndelVariantAnnotationsModel.calibration_set_scores]) }
    }
    if (defined(TrainIndelVariantAnnotationsModel.negative_model_scorer_pickle)) {
        call FF.FinalizeToFile as FinalizeIndelTrainVariantAnnotationsNegativeModelScorer { input: outdir = recalibration_model_dir, keyfile = keyfile, file = select_first([TrainIndelVariantAnnotationsModel.negative_model_scorer_pickle]) }
    }

    # ScoreVariantAnnotations
    # This was done per-contig, so we need to finalize per-contig:
    scatter (idx_3 in range(length(MakeChrIntervalList.contig_interval_list_files))) {

        String contig_3 = MakeChrIntervalList.chrs[idx_3][0]

        call FF.FinalizeToFile as FinalizeScoreSnpVariantAnnotationsScoredVcf { input: outdir = recalibration_results_dir + "/" + contig_3, keyfile = keyfile, file = ScoreSnpVariantAnnotations.scored_vcf[idx_3] }
        call FF.FinalizeToFile as FinalizeScoreSnpVariantAnnotationsScoredVcfIndex { input: outdir = recalibration_results_dir + "/" + contig_3, keyfile = keyfile, file = ScoreSnpVariantAnnotations.scored_vcf_index[idx_3] }
        if (defined(ScoreSnpVariantAnnotations.annotations_hdf5)) {
            call FF.FinalizeToFile as FinalizeScoreSnpVariantAnnotationsAnnotationsHdf5 { input: outdir = recalibration_results_dir + "/" + contig_3, keyfile = keyfile, file = select_first([ScoreSnpVariantAnnotations.annotations_hdf5[idx_3]]) }
        }
        if (defined(ScoreSnpVariantAnnotations.scores_hdf5)) {
            call FF.FinalizeToFile as FinalizeScoreSnpVariantAnnotationsScoresHdf5 { input: outdir = recalibration_results_dir + "/" + contig_3, keyfile = keyfile, file = select_first([ScoreSnpVariantAnnotations.scores_hdf5[idx_3]]) }
        }

        call FF.FinalizeToFile as FinalizeScoreIndelVariantAnnotationsScoredVcf { input: outdir = recalibration_results_dir + "/" + contig_3, keyfile = keyfile, file = ScoreIndelVariantAnnotations.scored_vcf[idx_3] }
        call FF.FinalizeToFile as FinalizeScoreIndelVariantAnnotationsScoredVcfIndex { input: outdir = recalibration_results_dir + "/" + contig_3, keyfile = keyfile, file = ScoreIndelVariantAnnotations.scored_vcf_index[idx_3] }
        if (defined(ScoreIndelVariantAnnotations.annotations_hdf5)) {
            call FF.FinalizeToFile as FinalizeScoreIndelVariantAnnotationsAnnotationsHdf5 { input: outdir = recalibration_results_dir + "/" + contig_3, keyfile = keyfile, file = select_first([ScoreIndelVariantAnnotations.annotations_hdf5[idx_3]]) }
        }
        if (defined(ScoreIndelVariantAnnotations.scores_hdf5)) {
            call FF.FinalizeToFile as FinalizeScoreIndelVariantAnnotationsScoresHdf5 { input: outdir = recalibration_results_dir + "/" + contig_3, keyfile = keyfile, file = select_first([ScoreIndelVariantAnnotations.scores_hdf5[idx_3]]) }
        }
    }

    # Make an alias for the functionally annotated data:
    if (defined(snpeff_db)) {
        File annotated_vcf = FinalizeVETSVCF.gcs_path
        File annotated_vcf_tbi = FinalizeVETSTBI.gcs_path
    }

    output {
        String genomicsDB = FinalizeGenomicsDB.gcs_dir

        File raw_joint_vcf     = FinalizeRawVCF.gcs_path
        File raw_joint_vcf_tbi = FinalizeRawTBI.gcs_path

        File joint_recalibrated_vcf     = FinalizeVETSVCF.gcs_path
        File joint_recalibrated_vcf_tbi = FinalizeVETSTBI.gcs_path

        File? annotated_joint_vcf     = annotated_vcf
        File? annotated_joint_vcf_tbi = annotated_vcf_tbi

        File joint_mt = CreateHailMatrixTable.gcs_path
        File joint_zarr = ConvertToZarr.gcs_path
    }
}


