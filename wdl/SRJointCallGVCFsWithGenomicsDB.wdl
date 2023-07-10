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

        Float snp_filter_level = 99.7
        Array[String] snp_recalibration_annotation_values = ["QD", "FS", "SOR", "MQRankSum", "ReadPosRankSum"]
        Array[Float] snp_recalibration_tranche_values = [100.0, 99.95, 99.9, 99.8, 99.6, 99.5, 99.4, 99.3, 99.0, 98.0, 97.0, 90.0 ]

        Array[File] snp_known_reference_variants
        Array[File] snp_known_reference_variants_index
        Array[File] snp_known_reference_variants_identifier
        Array[Boolean] snp_is_known
        Array[Boolean] snp_is_training
        Array[Boolean] snp_is_truth
        Array[Float] snp_prior
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
        Array[Float] indel_prior
        Int indel_max_gaussians = 8

        Array[File]?   annotation_bed_files
        Array[File]?   annotation_bed_file_indexes
        Array[String]? annotation_bed_file_annotation_names

        File? snpeff_db

        String prefix

        Boolean convert_to_zarr = false

        String gcs_out_root_dir
    }

    parameter_meta {
        gvcfs:            "GCS paths to gVCF files"
        gvcf_indices:     "GCS paths to gVCF tbi files"
        ref_map_file:     "table indicating reference sequence and auxillary file locations"
        prefix:           "prefix for output joint-called gVCF and tabix index"
        gcs_out_root_dir: "GCS bucket to store the reads, variants, and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/LRJointCallGVCFsWithGenomicsDB/~{prefix}"

    Map[String, String] ref_map = read_map(ref_map_file)

    # From WARP:
    # For small callsets (fewer than 1000 samples) we can gather the VCF shards and collect metrics directly.
    # For anything larger, we need to keep the VCF sharded and gather metrics collected from them.
    # We allow overriding this default behavior for testing / special requests.
    Boolean is_small_callset = length(gvcfs) <= 1000

    # Create interval list over which to shard the processing:
    call UTILS.MakeChrIntervalList as MakeChrIntervalList {
        input:
            ref_dict = ref_map['dict'],
    }

    # Reblock our GVCFs:
    scatter (idx_1 in range(length(gvcfs))) {
        call SRJOINT.ReblockGVCF as ReblockGVCFs {
            input:
                input_gvcf = gvcfs[idx_1],
                input_gvcf_index = gvcf_indices[idx_1],
                ref_fasta       = ref_map['fasta'],
                ref_fasta_fai   = ref_map['fai'],
                ref_dict        = ref_map['dict'],
                # Get rid of any and all suffixes:
                prefix = basename(basename(basename(gvcfs[idx_1], ".g.vcf.gz"), ".vcf.gz"), ".vcf")
        }
    }

    # Create sample-name map:
    call SRJOINT.CreateSampleNameMap as CreateSampleNameMap {
        input:
            gvcfs = gvcfs,
            prefix = prefix
    }

    # Shard by contig for speed:
    scatter (idx_2 in range(length(MakeChrIntervalList.contig_interval_list_files))) {

        String contig = MakeChrIntervalList.chrs[idx_2][0]
        File contig_interval_list = MakeChrIntervalList.contig_interval_list_files[idx_2]

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

    # Now we run VariantRecalibrator for indels and snps:
    call VARUTIL.IndelsVariantRecalibrator as TrainVQSROnHCIndelVariants {
        input:
            vcfs = MakeSitesOnlyVCF.sites_only_vcf,
            vcf_indices = MakeSitesOnlyVCF.sites_only_vcf_index,
            prefix = prefix + ".indels",
            recalibration_tranche_values = indel_recalibration_tranche_values,
            recalibration_annotation_values = indel_recalibration_annotation_values,
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
            vcfs = MakeSitesOnlyVCF.sites_only_vcf,
            vcf_indices = MakeSitesOnlyVCF.sites_only_vcf_index,
            prefix = prefix + ".snps",
            recalibration_tranche_values = snp_recalibration_tranche_values,
            recalibration_annotation_values = snp_recalibration_annotation_values,
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

    # Shard by contig for speed:
    scatter (idx_3 in range(length(JointCallGVCFs.output_vcf))) {

        File joint_called_vcf = JointCallGVCFs.output_vcf[idx_3]
        File joint_called_vcf_index = JointCallGVCFs.output_vcf[idx_3]

        call VARUTIL.ApplyVqsr as ApplyVqsr {
            input:
                vcf = joint_called_vcf,
                vcf_index = joint_called_vcf_index,

                prefix = basename(joint_called_vcf, ".vcf") + ".vqsr",

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
                    prefix = basename(ApplyVqsr.recalibrated_vcf, ".vcf") + ".region_annotated"
            }
        }

        File recalibrated_vcf = select_first([AnnotateVcfRegions.annotated_vcf, ApplyVqsr.recalibrated_vcf])
        File recalibrated_vcf_index = select_first([AnnotateVcfRegions.annotated_vcf_index, ApplyVqsr.recalibrated_vcf_index])

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
    call VARUTIL.GatherVcfs as GatherRecalibratedVcfs {
        input:
            input_vcfs = vcf_for_merging,
            input_vcf_indices = vcf_index_for_merging,
            prefix = prefix + ".recalibrated.combined"
    }

    # Convert to Zarr
    call SGKit.ConvertToZarrStore as ConvertToZarr {
        input:
            gvcf = GatherRecalibratedVcfs.output_vcf,
            tbi = GatherRecalibratedVcfs.output_vcf_index,
            prefix = prefix,
            outdir = outdir
    }

    # Convert the output to a HAIL Matrix Table:
    call Hail.ConvertToHailMT as CreateHailMatrixTable {
        input:
            gvcf = GatherRecalibratedVcfs.output_vcf,
            tbi = GatherRecalibratedVcfs.output_vcf_index,
            reference = sub(sub(ref_map["fasta"], "^.*/", ""), "\.[fasta]*$", ""),
            ref_fasta = ref_map["fasta"],
            ref_fai = ref_map["fai"],
            prefix = prefix,
            outdir = outdir
    }

    # Finalize:
    File keyfile = CreateHailMatrixTable.monitoring_log

#    call FF.FinalizeToDir as FinalizeGenomicsDB { input: outdir = outdir + "/GenomicsDB", keyfile = keyfile, file = ImportGVCFsIntoGenomicsDB.output_genomicsdb }

    call FF.FinalizeToFile as FinalizeRawVCF { input: outdir = outdir, keyfile = keyfile, file = GatherRawVcfs.output_vcf }
    call FF.FinalizeToFile as FinalizeRawTBI { input: outdir = outdir, keyfile = keyfile, file = GatherRawVcfs.output_vcf_index }

#    call FF.FinalizeToDir as FinalizeIndelRecalFile { input: outdir = outdir + "/recalibration_files", keyfile = keyfile, file = TrainVQSROnHCIndelVariants.recalibration }
#    call FF.FinalizeToDir as FinalizeIndelRecalIndex { input: outdir = outdir + "/recalibration_files, keyfile = keyfile, file = TrainVQSROnHCIndelVariants.recalibration_index }
#    call FF.FinalizeToDir as FinalizeIndelRecalTranches { input: outdir = outdir + "/recalibration_files, keyfile = keyfile, file = TrainVQSROnHCIndelVariants.tranches }
#    call FF.FinalizeToDir as FinalizeIndelRecalModelReport { input: outdir = outdir + "/recalibration_files, keyfile = keyfile, file = TrainVQSROnHCIndelVariants.model_report }

#    call FF.FinalizeToDir as FinalizeSnpRecalFile { input: outdir = outdir + "/recalibration_files, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.recalibration }
#    call FF.FinalizeToDir as FinalizeSnpRecalIndex { input: outdir = outdir + "/recalibration_files, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.recalibration_index }
#    call FF.FinalizeToDir as FinalizeSnpRecalTranches { input: outdir = outdir + "/recalibration_files, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.tranches }
#    call FF.FinalizeToDir as FinalizeSnpRecalModelReport { input: outdir = outdir + "/recalibration_files, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.model_report }

    call FF.FinalizeToFile as FinalizeVQSRVCF { input: outdir = outdir, keyfile = keyfile, file = GatherRecalibratedVcfs.output_vcf }
    call FF.FinalizeToFile as FinalizeVQSRTBI { input: outdir = outdir, keyfile = keyfile, file = GatherRecalibratedVcfs.output_vcf_index }

#    if (defined(annotation_bed_files)) {
#        call FF.FinalizeToFile as FinalizeRegionAnnotatedVcf { input: outdir = outdir, keyfile = keyfile, file = select_first([AnnotateVcfRegions.annotated_vcf]) }
#        call FF.FinalizeToFile as FinalizeRegionAnnotatedVcfIndex { input: outdir = outdir, keyfile = keyfile, file = select_first([AnnotateVcfRegions.annotated_vcf_index]) }
#    }

    ##########
    # store the results into designated bucket
    ##########

    output {
#        File genomicsDB = FinalizeGenomicsDB.gcs_path

        File raw_joint_vcf     = FinalizeRawVCF.gcs_path
        File raw_joint_vcf_tbi = FinalizeRawTBI.gcs_path

#        Array[File?] vqsr_indel_recal_file         = FinalizeIndelRecalFile.gcs_path
#        Array[File?] vqsr_indel_recal_file_index   = FinalizeIndelRecalIndex.gcs_path
#        Array[File?] vqsr_indel_recal_tranches     = FinalizeIndelRecalTranches.gcs_path
#        Array[File?] vqsr_indel_recal_model_report = FinalizeIndelRecalModelReport.gcs_path
#
#        Array[File?] vqsr_snp_recal_file         = FinalizeSnpRecalFile.gcs_path
#        Array[File?] vqsr_snp_recal_file_index   = FinalizeSnpRecalIndex.gcs_path
#        Array[File?] vqsr_snp_recal_tranches     = FinalizeSnpRecalTranches.gcs_path
#        Array[File?] vqsr_snp_recal_model_report = FinalizeSnpRecalModelReport.gcs_path

        File joint_recalibrated_vcf     = FinalizeVQSRVCF.gcs_path
        File joint_recalibrated_vcf_tbi = FinalizeVQSRTBI.gcs_path

#        File? annotated_joint_vcf     = AnnotateVcfRegions.annotated_vcf
#        File? annotated_joint_vcf_tbi = AnnotateVcfRegions.annotated_vcf_index

        File joint_mt = CreateHailMatrixTable.gcs_path
        File joint_zarr = ConvertToZarr.gcs_path
    }
}


