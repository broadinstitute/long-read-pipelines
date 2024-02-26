version 1.0

import "../../tasks/VariantCalling/SRJointGenotyping.wdl" as SRJOINT
import "../../tasks/Utility/VariantUtils.wdl" as VARUTIL
import "../../tasks/Utility/Utils.wdl" as UTILS
import "../../tasks/Utility/Finalize.wdl" as FF
import "../../tasks/Z_One_Off_Analyses/Pf_Niare_HaplotypeCaller.wdl" as Niare_HC

workflow SRJointCallGVCFsWithGenomicsDB_Pf_Niare_VQSR {

    meta {
        author: "Jonn Smith"
        description: "This workflow implements a modified version of the joint calling pipeline from Niare et al. (https://doi.org/10.1186/s12936-023-04632-0) using LRMA conventions.  The modification is that this pipeline uses VETS instead of VQSR."
    }
    parameter_meta {
        gvcfs:  "Array of GVCF files to use as inputs for joint calling."
        gvcf_indices:   "Array of gvcf index files for `gvcfs`.  Order should correspond to that in `gvcfs`."
        ref_map_file:  "Reference map file indicating reference sequence and auxillary file locations"

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

        prefix: "Prefix to use for output files."
        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs."
    }

    input {
        Array[File] gvcfs
        Array[File] gvcf_indices

        File ref_map_file

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

        String prefix

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/SRJointCallGVCFsWithGenomicsDB_Pf_Niare_VQSR/~{prefix}"

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

        call Niare_HC.GenomicsDbImport as GenomicsDbImport {
            input:
                sample_name_map = CreateSampleNameMap.sample_name_map,
                interval_list = contig_interval_list,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                prefix = prefix + "." + contig,
        }

        # Shard again by contig chunk:
        call UTILS.SplitContigToIntervals as SplitContigToIntervals {
            input:
                ref_dict = ref_map['dict'],
                contig = contig
        }

        scatter (idx_2 in range(length(SplitContigToIntervals.individual_bed_files))) {

            File genotype_gvcfs_intervals = SplitContigToIntervals.individual_bed_files[idx_2]

            call Niare_HC.GenotypeGVCFs as GenotypeGVCFs {
                input:
                    input_gvcf_data = GenomicsDbImport.output_genomicsdb,
                    interval_list = genotype_gvcfs_intervals,
                    ref_fasta         = ref_map['fasta'],
                    ref_fasta_fai     = ref_map['fai'],
                    ref_dict          = ref_map['dict'],
                    prefix = prefix + "." + contig + ".raw",
            }
        }

        # Merge all raw VCFs:
        call VARUTIL.GatherVcfs as GatherVcfs {
            input:
                input_vcfs = GenotypeGVCFs.output_vcf,
                input_vcf_indices = GenotypeGVCFs.output_vcf_index,
                prefix = prefix + "." + contig + ".raw.merged",
        }

        # First make a sites-only VCF for recal (smaller file, easier to work with):
        call VARUTIL.MakeSitesOnlyVcf as MakeSitesOnlyVCF {
            input:
                vcf = GatherVcfs.output_vcf,
                vcf_index = GatherVcfs.output_vcf_index,
                prefix = prefix + "." + contig + ".sites_only"
        }
    }

    # Merge all sites-only VCFs
    call VARUTIL.GatherVcfs as GatherSitesOnlyVCFs {
        input:
            input_vcfs = MakeSitesOnlyVCF.sites_only_vcf,
            input_vcf_indices = MakeSitesOnlyVCF.sites_only_vcf_index,
            prefix = prefix + ".sites_only"
    }

    ########################################################################
    # Call VETS / VQSR-lite:
    call VARUTIL.ExtractVariantAnnotations as ExtractIndelVariantAnnotations {
        input:
            vcf = GatherSitesOnlyVCFs.output_vcf,
            vcf_index = GatherSitesOnlyVCFs.output_vcf_index,

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
            vcf = GatherSitesOnlyVCFs.output_vcf,
            vcf_index = GatherSitesOnlyVCFs.output_vcf_index,

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
    scatter (idx_3 in range(length(GatherVcfs.output_vcf))) {

        String contig_2 = MakeChrIntervalList.chrs[idx_3][0]
        File joint_called_vcf = GatherVcfs.output_vcf[idx_3]
        File joint_called_vcf_index = GatherVcfs.output_vcf_index[idx_3]

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

        File vcf_for_merging = select_first([AnnotateVcfRegions.annotated_vcf, ScoreIndelVariantAnnotations.scored_vcf])
        File vcf_index_for_merging = select_first([AnnotateVcfRegions.annotated_vcf_index, ScoreIndelVariantAnnotations.scored_vcf_index])
    }

    # Consolidate files:
    call VARUTIL.GatherVcfs as GatherRescoredVcfs {
        input:
            input_vcfs = vcf_for_merging,
            input_vcf_indices = vcf_index_for_merging,
            prefix = prefix + ".rescored.combined"
    }

    ################################
    # Finalize the regular output files:
    ############
    File keyfile = GatherRescoredVcfs.output_vcf_index

    call FF.FinalizeToFile as FinalizeVETSVCF { input: outdir = outdir, keyfile = keyfile, file = GatherRescoredVcfs.output_vcf }
    call FF.FinalizeToFile as FinalizeVETSTBI { input: outdir = outdir, keyfile = keyfile, file = GatherRescoredVcfs.output_vcf_index }


    output {
        File joint_recalibrated_vcf     = FinalizeVETSVCF.gcs_path
        File joint_recalibrated_vcf_tbi = FinalizeVETSTBI.gcs_path
    }
}

