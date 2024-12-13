version 1.0

import "../../../tasks/VariantCalling/SRJointGenotyping.wdl" as SRJOINT
import "../../../tasks/Utility/VariantUtils.wdl" as VARUTIL
import "../../../tasks/Utility/Utils.wdl" as UTILS
import "../../../tasks/Utility/Hail.wdl" as Hail
import "../../../tasks/TertiaryAnalysis/FunctionalAnnotation.wdl" as FUNK
import "../../../tasks/Utility/SGKit.wdl" as SGKit
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow SRJointCallGVCFsWithGenomicsDB {

    meta {
        author: "Jonn Smith"
        description: "A workflow that performs joint calling on single-sample gVCFs from GATK4 HaplotypeCaller using GenomicsDB."
    }
    parameter_meta {
        gvcfs:  "Array of GVCF files to use as inputs for joint calling."
        gvcf_indices:   "Array of gvcf index files for `gvcfs`.  Order should correspond to that in `gvcfs`."
        ref_map_file:  "Reference map file indicating reference sequence and auxillary file locations"

        heterozygosity: "Joint Genotyping Parameter - Heterozygosity value used to compute prior likelihoods for any locus. See the GATKDocs for full details on the meaning of this population genetics concept"
        heterozygosity_stdev: "Joint Genotyping Parameter - Standard deviation of heterozygosity for SNP and indel calling."
        indel_heterozygosity: "Joint Genotyping Parameter - Heterozygosity for indel calling. See the GATKDocs for heterozygosity for full details on the meaning of this population genetics concept"

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

        annotation_bed_files:   "Array of bed files to use to FILTER/annotate variants in the output file.  Annotations will be placed in the FILTER column, effectively filtering variants that overlap these regions."
        annotation_bed_file_indexes:    "Array of bed indexes for `annotation_bed_files`.  Order should correspond to `annotation_bed_files`."
        annotation_bed_file_annotation_names:   "Array of names/FILTER column entries to use for each given file in `annotation_bed_files`.  Order should correspond to `annotation_bed_files`."

        use_gnarly_genotyper: "If true, the `GnarlyGenotyper` will be used, greatly speeding up joint genotyping (at the cost of potentially lower accuracy).  Setting this to true is recommended for large callsets.  If false, `GenotypeGVCFs` will be used to generate the final VCF.  Default is false."

        shard_max_interval_size_bp: "Maximum size of the interval on each shard.  This along with the given sequence dictionary determines how many shards there will be.  To shard by contig, set to a very high number.  Default is 999999999."

        prefix: "Prefix to use for output files."

        background_sample_gvcfs: "Array of GVCFs to use as background samples for joint calling."
        background_sample_gvcf_indices: "Array of GVCF index files for `background_sample_gvcfs`.  Order should correspond to that in `background_sample_gvcfs`."

        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs.  If no bucket is given, outputs will not be finalized and instead will remain in their native execution location."
    }

    input {
        Array[File] gvcfs
        Array[File] gvcf_indices

        File ref_map_file

        Float heterozygosity = 0.001
        Float heterozygosity_stdev = 0.01
        Float indel_heterozygosity = 0.000125

        Float snp_calibration_sensitivity = 0.99
        Int snp_max_unlabeled_variants = 0
        # TODO: Fix the annotations here to include the missing ones.  Must debug.
#        Array[String] snp_recalibration_annotation_values = [ "BaseQRankSum", "ExcessHet", "FS", "HAPCOMP", "HAPDOM", "HEC", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR", "DP" ]
        Array[String] snp_recalibration_annotation_values = [ "BaseQRankSum", "ExcessHet", "FS", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR", "DP" ]

        Array[File] snp_known_reference_variants
        Array[File] snp_known_reference_variants_index
        Array[File] snp_known_reference_variants_identifier
        Array[Boolean] snp_is_training
        Array[Boolean] snp_is_calibration

        Float indel_calibration_sensitivity = 0.99
        Int indel_max_unlabeled_variants = 0
        # TODO: Fix the annotations here to include the missing ones.  Must debug.
#        Array[String] indel_recalibration_annotation_values = [ "BaseQRankSum", "ExcessHet", "FS", "HAPCOMP", "HAPDOM", "HEC", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR", "DP" ]
        Array[String] indel_recalibration_annotation_values = [ "BaseQRankSum", "ExcessHet", "FS", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR", "DP" ]

        Array[File] indel_known_reference_variants
        Array[File] indel_known_reference_variants_index
        Array[File] indel_known_reference_variants_identifier
        Array[Boolean] indel_is_training
        Array[Boolean] indel_is_calibration

        Array[File]?   annotation_bed_files
        Array[File]?   annotation_bed_file_indexes
        Array[String]? annotation_bed_file_annotation_names

        File? snpeff_db
        String? snpeff_db_identifier

        File? interval_list

        Boolean use_gnarly_genotyper = false

        Int shard_max_interval_size_bp = 999999999

        String prefix

        Array[Array[File]]? background_sample_gvcfs
        Array[Array[File]]? background_sample_gvcf_indices

        String? gcs_out_root_dir
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    # Resolve the db_snp_vcf file, with preference to the db_snp_vcf file if it exists:
    call UTILS.ResolveMapKeysInPriorityOrder as ResolveMapKeysInPriorityOrder {
        input:
            map = ref_map,
            keys = ["test_bad_key_should_not_be_found", "dbsnp_vcf", "known_sites_vcf"]
    }
    File db_snp_vcf = ref_map[ResolveMapKeysInPriorityOrder.key]

    # Create sample-name map:
    call SRJOINT.CreateSampleNameMap as CreateSampleNameMap {
        input:
            gvcfs = gvcfs,
            background_sample_gvcfs = if defined(background_sample_gvcfs) then flatten(select_first([background_sample_gvcfs])) else [],
            prefix = prefix
    }

    # Get our interval list:
    if (!defined(interval_list)) {
        # If we have to, create interval list over which to shard the processing:
        call UTILS.MakeIntervalListFromSequenceDictionary as MakeIntervalListFromSequenceDictionary {
            input:
                ref_dict = ref_map['dict'],
                max_interval_size = shard_max_interval_size_bp
        }
    }
    File actual_interval_list = select_first([interval_list, MakeIntervalListFromSequenceDictionary.interval_list])

    # Get the interval name info for our files below:
    call UTILS.ExtractIntervalNamesFromIntervalOrBamFile as ExtractIntervalNamesFromIntervalOrBamFile {
        input:
            interval_file = actual_interval_list
    }

    # Shard by contig for speed:
    scatter (idx_1 in range(length(ExtractIntervalNamesFromIntervalOrBamFile.interval_info))) {

        String interval_name = ExtractIntervalNamesFromIntervalOrBamFile.interval_info[idx_1][0] + "_" + ExtractIntervalNamesFromIntervalOrBamFile.interval_info[idx_1][1] + "_" + ExtractIntervalNamesFromIntervalOrBamFile.interval_info[idx_1][2]

        # To make sure the interval names and the files themselves correspond, we need to make the
        # interval list file here:
        call UTILS.CreateIntervalListFileFromIntervalInfo as CreateIntervalListFileFromIntervalInfo {
            input:
                contig = ExtractIntervalNamesFromIntervalOrBamFile.interval_info[idx_1][0],
                start = ExtractIntervalNamesFromIntervalOrBamFile.interval_info[idx_1][1],
                end = ExtractIntervalNamesFromIntervalOrBamFile.interval_info[idx_1][2]
        }

        # Import our data into GenomicsDB:
        call SRJOINT.ImportGVCFs as ImportGVCFsIntoGenomicsDB {
            input:
                sample_name_map = CreateSampleNameMap.sample_name_map,
                interval_list   = CreateIntervalListFileFromIntervalInfo.interval_list,
                ref_fasta       = ref_map['fasta'],
                ref_fasta_fai   = ref_map['fai'],
                ref_dict        = ref_map['dict'],
                prefix          = prefix + "." + interval_name,
                batch_size      = 50,
                # We need to override this because we're not actually sending the GVCF over (just a list)
                # ALSO, we're currently tarring the genomicsDB, so we need at least double the space here, plus some slop:
                runtime_attr_override = object {disk_gb: 10 + (3 * CreateSampleNameMap.total_gvcf_size_gb) + (2 * ceil(size(ref_map['fasta'], "GB"))), preemptible_tries: 0}
        }

        # Joint call
        if (use_gnarly_genotyper) {
            call SRJOINT.GnarlyGenotypeGVCFs as GnarlyJointCallGVCFs {
                input:
                    input_gvcf_data = ImportGVCFsIntoGenomicsDB.output_genomicsdb,
                    interval_list   = CreateIntervalListFileFromIntervalInfo.interval_list,
                    ref_fasta       = ref_map['fasta'],
                    ref_fasta_fai   = ref_map['fai'],
                    ref_dict        = ref_map['dict'],
                    dbsnp_vcf       = db_snp_vcf,
                    prefix          = prefix + "." + interval_name + ".gnarly_genotyper.raw",
                    heterozygosity = heterozygosity,
                    heterozygosity_stdev = heterozygosity_stdev,
                    indel_heterozygosity = indel_heterozygosity,
                    runtime_attr_override = object {preemptible_tries: 0},  # Disable preemption for prototype.
            }
        }
        if (!use_gnarly_genotyper) {
            call SRJOINT.GenotypeGVCFs as JointCallGVCFs {
                input:
                    input_gvcf_data = ImportGVCFsIntoGenomicsDB.output_genomicsdb,
                    interval_list   = CreateIntervalListFileFromIntervalInfo.interval_list,
                    ref_fasta       = ref_map['fasta'],
                    ref_fasta_fai   = ref_map['fai'],
                    ref_dict        = ref_map['dict'],
                    dbsnp_vcf       = db_snp_vcf,
                    prefix          = prefix + "." + interval_name + ".genotype_gvcfs.raw",
                    heterozygosity = heterozygosity,
                    heterozygosity_stdev = heterozygosity_stdev,
                    indel_heterozygosity = indel_heterozygosity,
                    runtime_attr_override = object {preemptible_tries: 0},  # Disable preemption for prototype.
            }
        }
        # Select the VCF + index for the raw joint called file:
        File joint_vcf = select_first([GnarlyJointCallGVCFs.output_vcf, JointCallGVCFs.output_vcf])
        File joint_vcf_index = select_first([GnarlyJointCallGVCFs.output_vcf_index, JointCallGVCFs.output_vcf_index])

        # First make a sites-only VCF for recal (smaller file, easier to work with):
        call VARUTIL.MakeSitesOnlyVcf as MakeSitesOnlyVCF {
            input:
                vcf = joint_vcf,
                vcf_index = joint_vcf_index,
                prefix = prefix + "." + interval_name + ".sites_only"
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
            annotation_hdf5 = ExtractSnpVariantAnnotations.annotation_hdf5,
            mode = "SNP",
            prefix = prefix,
    }

    # Shard by contig for speed:
    scatter (idx_2 in range(length(ExtractIntervalNamesFromIntervalOrBamFile.interval_info))) {

        String interval_name2 = ExtractIntervalNamesFromIntervalOrBamFile.interval_info[idx_2][0] + "_" + ExtractIntervalNamesFromIntervalOrBamFile.interval_info[idx_2][1] + "_" + ExtractIntervalNamesFromIntervalOrBamFile.interval_info[idx_2][2]
        File joint_called_vcf = joint_vcf[idx_2]
        File joint_called_vcf_index = joint_vcf_index[idx_2]

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
                prefix = prefix + "_SNP_" + interval_name2,
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
                prefix = prefix + "_ALL_" + interval_name2,
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
                    snpeff_db = select_first([snpeff_db]),
                    snpeff_db_identifier = select_first([snpeff_db_identifier])
            }
        }

        File vcf_for_merging = select_first([FunctionallyAnnotate.annotated_vcf, recalibrated_vcf])
        File vcf_index_for_merging = select_first([FunctionallyAnnotate.annotated_vcf_index, recalibrated_vcf_index])
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
            vcf = GatherRescoredVcfs.output_vcf,
            tbi = GatherRescoredVcfs.output_vcf_index,
            prefix = prefix
    }

    # Convert the output to a HAIL Matrix Table:
    call Hail.ConvertToHailMT as CreateHailMatrixTable {
        input:
            gvcf = GatherRescoredVcfs.output_vcf,
            tbi = GatherRescoredVcfs.output_vcf_index,
            reference = sub(sub(ref_map["fasta"], "^.*/", ""), "\.[fasta]*$", ""),
            ref_fasta = ref_map["fasta"],
            ref_fai = ref_map["fai"],
            prefix = prefix
    }

    ################################
    # Finalize the regular output files:
    ############

    if (defined(gcs_out_root_dir)) {

        String concrete_gcs_out_root_dir = select_first([gcs_out_root_dir])
        String outdir = sub(concrete_gcs_out_root_dir, "/$", "") + "/SRJointCallGVCFsWithGenomicsDB/~{prefix}"

        String recalibration_dir = outdir + "/recalibration_files"
        String recalibration_model_dir = outdir + "/recalibration_files/model"
        String recalibration_results_dir = outdir + "/recalibration_files/results"
        String snpeff_results_dir = outdir + "/snpEff_results"

        File keyfile = CreateHailMatrixTable.completion_file

        call FF.FinalizeToDir as FinalizeGenomicsDB { input: outdir = outdir + "/GenomicsDB", keyfile = keyfile, files = ImportGVCFsIntoGenomicsDB.output_genomicsdb }

        call FF.FinalizeToFile as FinalizeVETSVCF { input: outdir = outdir, keyfile = keyfile, file = GatherRescoredVcfs.output_vcf }
        call FF.FinalizeToFile as FinalizeVETSTBI { input: outdir = outdir, keyfile = keyfile, file = GatherRescoredVcfs.output_vcf_index }

        if (defined(snpeff_db)) {
            call FF.FinalizeToDir as FinalizeSnpEffSummary { input: outdir = snpeff_results_dir, keyfile = keyfile, files = select_all(FunctionallyAnnotate.snpEff_summary) }
            call FF.FinalizeToDir as FinalizeSnpEffGenes { input: outdir = snpeff_results_dir, keyfile = keyfile, files = select_all(FunctionallyAnnotate.snpEff_genes) }
        }

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
        scatter (idx_3 in range(length(ExtractIntervalNamesFromIntervalOrBamFile.interval_info))) {

            String interval_name3 = ExtractIntervalNamesFromIntervalOrBamFile.interval_info[idx_3][0] + "_" + ExtractIntervalNamesFromIntervalOrBamFile.interval_info[idx_3][1] + "_" + ExtractIntervalNamesFromIntervalOrBamFile.interval_info[idx_3][2]

            call FF.FinalizeToFile as FinalizeScoreSnpVariantAnnotationsScoredVcf { input: outdir = recalibration_results_dir + "/" + interval_name3, keyfile = keyfile, file = ScoreSnpVariantAnnotations.scored_vcf[idx_3] }
            call FF.FinalizeToFile as FinalizeScoreSnpVariantAnnotationsScoredVcfIndex { input: outdir = recalibration_results_dir + "/" + interval_name3, keyfile = keyfile, file = ScoreSnpVariantAnnotations.scored_vcf_index[idx_3] }
            if (defined(ScoreSnpVariantAnnotations.annotations_hdf5)) {
                call FF.FinalizeToFile as FinalizeScoreSnpVariantAnnotationsAnnotationsHdf5 { input: outdir = recalibration_results_dir + "/" + interval_name3, keyfile = keyfile, file = select_first([ScoreSnpVariantAnnotations.annotations_hdf5[idx_3]]) }
            }
            if (defined(ScoreSnpVariantAnnotations.scores_hdf5)) {
                call FF.FinalizeToFile as FinalizeScoreSnpVariantAnnotationsScoresHdf5 { input: outdir = recalibration_results_dir + "/" + interval_name3, keyfile = keyfile, file = select_first([ScoreSnpVariantAnnotations.scores_hdf5[idx_3]]) }
            }

            call FF.FinalizeToFile as FinalizeScoreIndelVariantAnnotationsScoredVcf { input: outdir = recalibration_results_dir + "/" + interval_name3, keyfile = keyfile, file = ScoreIndelVariantAnnotations.scored_vcf[idx_3] }
            call FF.FinalizeToFile as FinalizeScoreIndelVariantAnnotationsScoredVcfIndex { input: outdir = recalibration_results_dir + "/" + interval_name3, keyfile = keyfile, file = ScoreIndelVariantAnnotations.scored_vcf_index[idx_3] }
            if (defined(ScoreIndelVariantAnnotations.annotations_hdf5)) {
                call FF.FinalizeToFile as FinalizeScoreIndelVariantAnnotationsAnnotationsHdf5 { input: outdir = recalibration_results_dir + "/" + interval_name3, keyfile = keyfile, file = select_first([ScoreIndelVariantAnnotations.annotations_hdf5[idx_3]]) }
            }
            if (defined(ScoreIndelVariantAnnotations.scores_hdf5)) {
                call FF.FinalizeToFile as FinalizeScoreIndelVariantAnnotationsScoresHdf5 { input: outdir = recalibration_results_dir + "/" + interval_name3, keyfile = keyfile, file = select_first([ScoreIndelVariantAnnotations.scores_hdf5[idx_3]]) }
            }
        }

        call FF.FinalizeToFile as FinalizeZarr { input: outdir = outdir, keyfile = keyfile, file = ConvertToZarr.zarr }

        call FF.FinalizeToFile as FinalizeHailMatrixTable {
            input:
                outdir = outdir,
                keyfile = keyfile,
                file = CreateHailMatrixTable.mt_tar
        }

        # Set up variable for outputs:
        Array[String] final_genomicsdb_location = [FinalizeGenomicsDB.gcs_dir]
    }

    # Make an alias for the functionally annotated data:
    if (defined(snpeff_db)) {
        File annotated_vcf = if defined(gcs_out_root_dir) then select_first([FinalizeVETSVCF.gcs_path]) else GatherRescoredVcfs.output_vcf
        File annotated_vcf_tbi = if defined(gcs_out_root_dir) then select_first([FinalizeVETSTBI.gcs_path]) else GatherRescoredVcfs.output_vcf_index

        Array[String] final_snpeff_summary = if defined(gcs_out_root_dir) then [select_first([FinalizeSnpEffSummary.gcs_dir])] else select_all(FunctionallyAnnotate.snpEff_summary)
        Array[String] final_snpEff_genes = if defined(gcs_out_root_dir) then [select_first([FinalizeSnpEffGenes.gcs_dir])] else select_all(FunctionallyAnnotate.snpEff_genes)
    }

    output {
        File joint_vcf_out     = select_first([FinalizeVETSVCF.gcs_path, GatherRescoredVcfs.output_vcf])
        File joint_vcf_out_tbi = select_first([FinalizeVETSTBI.gcs_path, GatherRescoredVcfs.output_vcf_index])

        File joint_mt = select_first([FinalizeHailMatrixTable.gcs_path, CreateHailMatrixTable.mt_tar])
        File joint_zarr = select_first([FinalizeZarr.gcs_path, ConvertToZarr.zarr])

        Array[String] genomicsDB = select_first([final_genomicsdb_location, ImportGVCFsIntoGenomicsDB.output_genomicsdb])

        Array[String]? snpEff_summary = final_snpeff_summary
        Array[String]? snpEff_genes = final_snpEff_genes
    }
}


