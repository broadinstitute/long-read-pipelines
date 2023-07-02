version 1.0

#############################################################################################################
## A workflow that performs joint calling on single-sample gVCFs from GATK4 HaplotypeCaller using GenomicsDB.
#############################################################################################################

import "tasks/SRJointGenotyping.wdl" as SRJOINT
import "tasks/VariantUtils.wdl" as VARUTIL
import "tasks/Finalize.wdl" as FF

workflow SRJointCallGVCFsWithGenomicsDB {
    input {
        Array[File] gvcfs
        Array[File] gvcf_indices

        File ref_map_file

        File interval_list

        Float snp_filter_level = 99.7
        Array[String] snp_recalibration_annotation_values = ["QD", "FS", "SOR", "MQRankSum", "ReadPosRankSum"]
        Array[Float] snp_recalibration_tranche_values = [100.0, 99.95, 99.9, 99.8, 99.6, 99.5, 99.4, 99.3, 99.0, 98.0, 97.0, 90.0 ]

        Array[File] snp_known_reference_variants
        Array[File] snp_known_reference_variants_index
        Array[File] snp_known_reference_variants_identifier
        Array[Boolean] snp_is_known
        Array[Boolean] snp_is_training
        Array[Boolean] snp_is_truth
        Array[Int] snp_prior
        Int snp_max_gaussians = 8

        Float indel_filter_level = 99.0
        Array[String] indel_recalibration_annotation_values = ["QD", "FS", "SOR", "MQRankSum", "ReadPosRankSum"]
        Array[Float] indel_recalibration_tranche_values = [100.0, 99.95, 99.9, 99.5, 99.0, 97.0, 96.0, 95.0, 94.0, 93.5, 93.0, 92.0, 91.0, 90.0]

        Array[File] indel_known_reference_variants
        Array[File] indel_known_reference_variants_index
        Array[File] indel_known_reference_variants_identifier
        Array[Boolean] indel_is_known
        Array[Boolean] indel_is_training
        Array[Boolean] indel_is_truth
        Array[Int] indel_prior
        Int indel_max_gaussians = 8

        Array[File]?   annotation_bed_files
        Array[File]?   annotation_bed_file_indexes
        Array[String]? annotation_bed_file_annotation_names

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

    # From WARP:
    # For small callsets (fewer than 1000 samples) we can gather the VCF shards and collect metrics directly.
    # For anything larger, we need to keep the VCF sharded and gather metrics collected from them.
    # We allow overriding this default behavior for testing / special requests.
    Boolean is_small_callset = length(gvcfs) <= 1000

    # Create sample-name map:
    call SRJOINT.CreateSampleNameMap as CreateSampleNameMap {
        input:
            gvcfs = gvcfs,
            prefix = prefix
    }

    # Import our data into GenomicsDB:
    call SRJOINT.ImportGVCFs as ImportGVCFsIntoGenomicsDB {
        input:
            sample_name_map = CreateSampleNameMap.sample_name_map,
            interval_list   = interval_list,
            ref_fasta       = ref_map['fasta'],
            ref_fasta_fai   = ref_map['fai'],
            ref_dict        = ref_map['dict'],
            prefix          = prefix,
            batch_size      = 50,
    }

    # Joint call
    call SRJOINT.GenotypeGVCFs as JointCallGVCFs {
        input:
            input_gvcf_data = ImportGVCFsIntoGenomicsDB.output_genomicsdb,
            interval_list   = interval_list,
            ref_fasta       = ref_map['fasta'],
            ref_fasta_fai   = ref_map['fai'],
            ref_dict        = ref_map['dict'],
            dbsnp_vcf       = ref_map["known_sites_vcf"],
            prefix          = prefix,
    }

    # First make a sites-only VCF for recal (smaller file, easier to work with):
    call VARUTIL.MakeSitesOnlyVcf as MakeSitesOnlyGVCF {
        input:
            vcf = JointCallGVCFs.output_vcf,
            vcf_index = JointCallGVCFs.output_vcf_index,
            prefix = prefix
    }

    # Now we run VariantRecalibrator for indels and snps:
    call VARUTIL.IndelsVariantRecalibrator as TrainVQSROnHCIndelVariants {
        input:
            vcf = MakeSitesOnlyGVCF.sites_only_vcf,
            vcf_index = MakeSitesOnlyGVCF.sites_only_vcf_index,
            prefix = prefix + ".indels",
            recalibration_tranche_values = indel_recalibration_tranche_values,
            recalibration_annotation_values = indel_recalibration_annotation_values,
#                known_reference_variants = [ref_map["known_sites_vcf"]],
#                known_reference_variants_index = [ref_map["known_sites_index"]],
#                known_reference_variants_identifier = ["pfcrosses"],
#                is_known = [true],
#                is_training = [true],
#                is_truth = [true],
#                prior = [15],
            known_reference_variants = indel_known_reference_variants,
            known_reference_variants_index = indel_known_reference_variants_index,
            known_reference_variants_identifier = indel_known_reference_variants_identifier,
            is_known = indel_is_known,
            is_training = indel_is_training,
            is_truth = indel_is_truth,
            prior = indel_prior,
            use_allele_specific_annotations = false,
            max_gaussians = indel_max_gaussians,
    }

    call VARUTIL.SNPsVariantRecalibratorCreateModel as TrainVQSROnHCSnpVariants {
        input:
            vcf = MakeSitesOnlyGVCF.sites_only_vcf,
            vcf_index = MakeSitesOnlyGVCF.sites_only_vcf_index,
            prefix = prefix + ".snps",
            recalibration_tranche_values = snp_recalibration_tranche_values,
            recalibration_annotation_values = snp_recalibration_annotation_values,
#                known_reference_variants = [ref_map["known_sites_vcf"]],
#                known_reference_variants_index = [ref_map["known_sites_index"]],
#                known_reference_variants_identifier = ["pfcrosses"],
#                is_known = [true],
#                is_training = [true],
#                is_truth = [true],
#                prior = [15],
            known_reference_variants = snp_known_reference_variants,
            known_reference_variants_index = snp_known_reference_variants_index,
            known_reference_variants_identifier = snp_known_reference_variants_identifier,
            is_known = snp_is_known,
            is_training = snp_is_training,
            is_truth = snp_is_truth,
            prior = snp_prior,
            use_allele_specific_annotations = false,
            max_gaussians = snp_max_gaussians,
    }

    call VARUTIL.ApplyVqsr as ApplyVqsr {
        input:
            vcf = JointCallGVCFs.output_vcf,
            vcf_index = JointCallGVCFs.output_vcf_index,

            prefix = prefix + ".vqsr_filtered",

            snps_recalibration = TrainVQSROnHCSnpVariants.recalibration,
            snps_recalibration_index = TrainVQSROnHCSnpVariants.recalibration_index,
            snps_tranches = TrainVQSROnHCSnpVariants.tranches,
            snp_filter_level = snp_filter_level,

            indels_recalibration = TrainVQSROnHCIndelVariants.recalibration,
            indels_recalibration_index = TrainVQSROnHCIndelVariants.recalibration_index,
            indels_tranches = TrainVQSROnHCIndelVariants.tranches,
            indel_filter_level = indel_filter_level,

            use_allele_specific_annotations = false,
    }

    # Now we need to annotate our variants by region:
    if (defined(annotation_bed_files)) {
        call VARUTIL.AnnotateVcfWithBedRegions as AnnotateVcfRegions {
            input:
                vcf = ApplyVqsr.recalibrated_vcf,
                vcf_index = ApplyVqsr.recalibrated_vcf_index,
                bed_files = select_first([annotation_bed_files]),
                bed_file_indexes = select_first([annotation_bed_file_indexes]),
                bed_file_annotation_names = select_first([annotation_bed_file_annotation_names]),
                prefix = prefix + ".region_annotated"
        }
    }

    # Finalize:
    File keyfile = select_first([AnnotateVcfRegions.annotated_vcf_index, ApplyVqsr.recalibrated_vcf_index])

    call FF.FinalizeToFile as FinalizeGenomicsDB { input: outdir = outdir, keyfile = keyfile, file = ImportGVCFsIntoGenomicsDB.output_genomicsdb }

    call FF.FinalizeToFile as FinalizeRawVCF { input: outdir = outdir, keyfile = keyfile, file = JointCallGVCFs.output_vcf }
    call FF.FinalizeToFile as FinalizeRawTBI { input: outdir = outdir, keyfile = keyfile, file = JointCallGVCFs.output_vcf_index }

    call FF.FinalizeToFile as FinalizeIndelRecalFile { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCIndelVariants.recalibration }
    call FF.FinalizeToFile as FinalizeIndelRecalIndex { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCIndelVariants.recalibration_index }
    call FF.FinalizeToFile as FinalizeIndelRecalTranches { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCIndelVariants.tranches }
    call FF.FinalizeToFile as FinalizeIndelRecalModelReport { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCIndelVariants.model_report }

    call FF.FinalizeToFile as FinalizeSnpRecalFile { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.recalibration }
    call FF.FinalizeToFile as FinalizeSnpRecalIndex { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.recalibration_index }
    call FF.FinalizeToFile as FinalizeSnpRecalTranches { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.tranches }
    call FF.FinalizeToFile as FinalizeSnpRecalModelReport { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.model_report }

    call FF.FinalizeToFile as FinalizeVQSRVCF { input: outdir = outdir, keyfile = keyfile, file = ApplyVqsr.recalibrated_vcf }
    call FF.FinalizeToFile as FinalizeVQSRTBI { input: outdir = outdir, keyfile = keyfile, file = ApplyVqsr.recalibrated_vcf_index }

    if (defined(annotation_bed_files)) {
        call FF.FinalizeToFile as FinalizeRegionAnnotatedVcf { input: outdir = outdir, keyfile = keyfile, file = select_first([AnnotateVcfRegions.annotated_vcf]) }
        call FF.FinalizeToFile as FinalizeRegionAnnotatedVcfIndex { input: outdir = outdir, keyfile = keyfile, file = select_first([AnnotateVcfRegions.annotated_vcf_index]) }
    }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File genomicsDB = FinalizeGenomicsDB.gcs_path

        File raw_joint_vcf     = FinalizeRawVCF.gcs_path
        File raw_joint_vcf_tbi = FinalizeRawTBI.gcs_path

        File? vqsr_indel_recal_file         = FinalizeIndelRecalFile.gcs_path
        File? vqsr_indel_recal_file_index   = FinalizeIndelRecalIndex.gcs_path
        File? vqsr_indel_recal_tranches     = FinalizeIndelRecalTranches.gcs_path
        File? vqsr_indel_recal_model_report = FinalizeIndelRecalModelReport.gcs_path

        File? vqsr_snp_recal_file         = FinalizeSnpRecalFile.gcs_path
        File? vqsr_snp_recal_file_index   = FinalizeSnpRecalIndex.gcs_path
        File? vqsr_snp_recal_tranches     = FinalizeSnpRecalTranches.gcs_path
        File? vqsr_snp_recal_model_report = FinalizeSnpRecalModelReport.gcs_path

        File joint_recalibrated_vcf     = FinalizeVQSRVCF.gcs_path
        File joint_recalibrated_vcf_tbi = FinalizeVQSRTBI.gcs_path
    }
}


