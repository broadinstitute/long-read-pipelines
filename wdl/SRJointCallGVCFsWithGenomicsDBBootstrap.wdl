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
    call UTILS.MakeChrIntervalList as t_001_MakeChrIntervalList {
        input:
            ref_dict = ref_map['dict'],
    }

    # Create sample-name map:
    call SRJOINT.CreateSampleNameMap as t_002_CreateSampleNameMap {
        input:
            gvcfs = gvcfs,
            prefix = prefix
    }

    # Shard by contig for speed:
    scatter (idx_1 in range(length(t_001_MakeChrIntervalList.contig_interval_list_files))) {

        String contig = t_001_MakeChrIntervalList.chrs[idx_1][0]
        File contig_interval_list = t_001_MakeChrIntervalList.contig_interval_list_files[idx_1]

        # Import our data into GenomicsDB:
        call SRJOINT.ImportGVCFs as t_003_ImportGVCFsIntoGenomicsDB {
            input:
                sample_name_map = t_002_CreateSampleNameMap.sample_name_map,
                interval_list   = contig_interval_list,
                ref_fasta       = ref_map['fasta'],
                ref_fasta_fai   = ref_map['fai'],
                ref_dict        = ref_map['dict'],
                prefix          = prefix + "." + contig,
                batch_size      = 50,
                # We need to override this because we're not actually sending the GVCF over (just a list)
                # ALSO, we're currently tarring the genomicsDB, so we need at least double the space here, plus some slop:
                runtime_attr_override = object {disk_gb: 10 + (3 * t_002_CreateSampleNameMap.total_gvcf_size_gb) + (2 * ceil(size(ref_map['fasta'], "GB"))), preemptible_tries: 0}
        }

        # Joint call
        call SRJOINT.GenotypeGVCFs as t_004_JointCallGVCFs {
            input:
                input_gvcf_data = t_003_ImportGVCFsIntoGenomicsDB.output_genomicsdb,
                interval_list   = contig_interval_list,
                ref_fasta       = ref_map['fasta'],
                ref_fasta_fai   = ref_map['fai'],
                ref_dict        = ref_map['dict'],
                # dbsnp_vcf       = ref_map["known_sites_vcf"], #Turn off dbsnp_vcf here
                prefix          = prefix + "." + contig + ".raw",
                runtime_attr_override = object {preemptible_tries: 0},  # Disable preemption for prototype.
        }

        # First make a sites-only VCF for recal (smaller file, easier to work with):
        call VARUTIL.MakeSitesOnlyVcf as t_005_MakeSitesOnlyVCF {
            input:
                vcf = t_004_JointCallGVCFs.output_vcf,
                vcf_index = t_004_JointCallGVCFs.output_vcf_index,
                prefix = prefix + "." + contig + ".sites_only"
        }
    }

    # Merge all sites-only VCFs
    call VARUTIL.GatherVcfs as t_006_MergeSitesOnlyVCFs {
        input:
            input_vcfs = t_005_MakeSitesOnlyVCF.sites_only_vcf,
            input_vcf_indices = t_005_MakeSitesOnlyVCF.sites_only_vcf_index,
            prefix = prefix + ".sites_only"
    }

    # Call hard filtering

    # Hard filter SNPs
    call VARUTIL.HardFilterVcfByGATKDefault_Snp as t_007_HardFilterSnps {
        input:
            vcf = t_006_MergeSitesOnlyVCFs.output_vcf,
            vcf_index = t_006_MergeSitesOnlyVCFs.output_vcf_index,
            ref_fasta = ref_map['fasta'],
            prefix = prefix
    }

    # Hard filter indels
    call VARUTIL.HardFilterVcfByGATKDefault_Indel as t_008_HardFilterIndels {
        input:
            vcf = t_006_MergeSitesOnlyVCFs.output_vcf,
            vcf_index = t_006_MergeSitesOnlyVCFs.output_vcf_index,
            ref_fasta = ref_map['fasta'],
            prefix = prefix
    }

    Array[File] filtered_vcfs = [t_007_HardFilterSnps.snp_filtered_vcf, t_008_HardFilterIndels.indel_filtered_vcf]

    # Merge Snp and Indel VCF files
    call VARUTIL.MergeAndSortVCFs as t_009_MergeSnpsAndIndels {
        input:
            vcfs = filtered_vcfs,
            ref_fasta_fai = ref_map["fai"],
            prefix = prefix,
            optional_flags = "--allow-overlaps"
    }

    # Shard by contig for speed:
    scatter (idx_2 in range(length(t_004_JointCallGVCFs.output_vcf))) {

        String contig_2 = t_001_MakeChrIntervalList.chrs[idx_2][0]
        File joint_called_vcf = t_004_JointCallGVCFs.output_vcf[idx_2]
        File joint_called_vcf_index = t_004_JointCallGVCFs.output_vcf_index[idx_2]

        # Annotating our variants by region:
        if (defined(annotation_bed_files)) {
            call VARUTIL.AnnotateVcfWithBedRegions as t_010_AnnotateVcfRegions {
                input:
                    vcf = t_009_MergeSnpsAndIndels.vcf,
                    vcf_index = t_009_MergeSnpsAndIndels.tbi,
                    bed_files = select_first([annotation_bed_files]),
                    bed_file_indexes = select_first([annotation_bed_file_indexes]),
                    bed_file_annotation_names = select_first([annotation_bed_file_annotation_names]),
                    prefix = basename(basename(t_009_MergeSnpsAndIndels.vcf, ".vcf.gz"), ".vcf") + ".region_annotated",
            }
        }

        File recalibrated_vcf = select_first([t_010_AnnotateVcfRegions.annotated_vcf, t_009_MergeSnpsAndIndels.vcf])
        File recalibrated_vcf_index = select_first([t_010_AnnotateVcfRegions.annotated_vcf_index, t_009_MergeSnpsAndIndels.tbi])

        # Now functionally annotate each VCF:
        if (defined(snpeff_db)) {
            call FUNK.FunctionallyAnnotateVariants as t_011_FunctionallyAnnotate {
                input:
                    vcf = recalibrated_vcf,
                    snpeff_db = select_first([snpeff_db])
            }
        }

        File vcf_for_merging = select_first([t_011_FunctionallyAnnotate.annotated_vcf, recalibrated_vcf])
        File vcf_index_for_merging = select_first([t_011_FunctionallyAnnotate.annotated_vcf_index, recalibrated_vcf_index])
    }

    # Consolidate files:
    call VARUTIL.GatherVcfs as t_012_GatherRawVcfs {
        input:
            input_vcfs = t_004_JointCallGVCFs.output_vcf,
            input_vcf_indices = t_004_JointCallGVCFs.output_vcf_index,
            prefix = prefix + ".raw.combined"
    }

    # Consolidate files:
    call VARUTIL.GatherVcfs as t_013_GatherFilteredVcfs {
        input:
            input_vcfs = vcf_for_merging,
            input_vcf_indices = vcf_index_for_merging,
            prefix = prefix + ".filtered.combined"
    }

    # Convert to Zarr
    call SGKit.ConvertToZarrStore as t_014_ConvertToZarr {
        input:
            gvcf = t_013_GatherFilteredVcfs.output_vcf,
            tbi = t_013_GatherFilteredVcfs.output_vcf_index,
            prefix = prefix,
            outdir = outdir
    }

    # Convert the output to a HAIL Matrix Table:
    call Hail.ConvertToHailMT as t_015_CreateHailMatrixTable {
        input:
            gvcf = t_013_GatherFilteredVcfs.output_vcf,
            tbi = t_013_GatherFilteredVcfs.output_vcf_index,
            reference = sub(sub(ref_map["fasta"], "^.*/", ""), "\.[fasta]*$", ""),
            ref_fasta = ref_map["fasta"],
            ref_fai = ref_map["fai"],
            prefix = prefix,
            outdir = outdir
    }

    ################################
    # Finalize the regular output files:
    ############

    File keyfile = t_015_CreateHailMatrixTable.completion_file
    String recalibration_dir = outdir + "/recalibration_files"
    String recalibration_model_dir = outdir + "/recalibration_files/model"
    String recalibration_results_dir = outdir + "/recalibration_files/results"
    String snpeff_results_dir = outdir + "/snpEff_results"

    call FF.FinalizeToDir as t_016_FinalizeGenomicsDB { input: outdir = outdir + "/GenomicsDB", keyfile = keyfile, files = t_003_ImportGVCFsIntoGenomicsDB.output_genomicsdb }

    call FF.FinalizeToFile as t_017_FinalizeRawVCF { input: outdir = outdir, keyfile = keyfile, file = t_012_GatherRawVcfs.output_vcf }
    call FF.FinalizeToFile as t_018_FinalizeRawTBI { input: outdir = outdir, keyfile = keyfile, file = t_012_GatherRawVcfs.output_vcf_index }

    call FF.FinalizeToFile as t_019_FinalizeFilteredVCF { input: outdir = outdir, keyfile = keyfile, file = t_013_GatherFilteredVcfs.output_vcf }
    call FF.FinalizeToFile as t_020_FinalizeFilteredTBI { input: outdir = outdir, keyfile = keyfile, file = t_013_GatherFilteredVcfs.output_vcf_index }

    if (defined(snpeff_db)) {
        call FF.FinalizeToDir as t_021_FinalizeSnpEffSummary { input: outdir = snpeff_results_dir, keyfile = keyfile, files = select_all(t_011_FunctionallyAnnotate.snpEff_summary) }
        call FF.FinalizeToDir as t_022_FinalizeSnpEffGenes { input: outdir = snpeff_results_dir, keyfile = keyfile, files = select_all(t_011_FunctionallyAnnotate.snpEff_genes) }
    }

    #####################################
    # Finalize the hard filtered files:
    #####################################
    call FF.FinalizeToFile as t_023_FinalizeFilteredSnpsVCF { input: outdir = outdir, keyfile = keyfile, file = t_007_HardFilterSnps.snp_filtered_vcf }
    call FF.FinalizeToFile as t_024_FinalizeFilteredSnpsTBI { input: outdir = outdir, keyfile = keyfile, file = t_007_HardFilterSnps.snp_filtered_vcf_index }

    call FF.FinalizeToFile as t_025_FinalizeFilteredIndelsVCF { input: outdir = outdir, keyfile = keyfile, file = t_008_HardFilterIndels.indel_filtered_vcf }
    call FF.FinalizeToFile as t_026_FinalizeFilteredIndelsTBI { input: outdir = outdir, keyfile = keyfile, file = t_008_HardFilterIndels.indel_filtered_vcf_index}


    # Make an alias for the functionally annotated data:
    if (defined(snpeff_db)) {
        File annotated_vcf = t_019_FinalizeFilteredVCF.gcs_path
        File annotated_vcf_tbi = t_020_FinalizeFilteredTBI.gcs_path
    }

    output {
        String genomicsDB = t_016_FinalizeGenomicsDB.gcs_dir

        File raw_joint_vcf     = t_017_FinalizeRawVCF.gcs_path
        File raw_joint_vcf_tbi = t_018_FinalizeRawTBI.gcs_path

        File joint_filtered_vcf     = t_019_FinalizeFilteredVCF.gcs_path
        File joint_filtered_vcf_tbi = t_019_FinalizeFilteredVCF.gcs_path

        File joint_mt = t_015_CreateHailMatrixTable.gcs_path
        File joint_zarr = t_014_ConvertToZarr.gcs_path

        File joint_snps_filtered_vcf     = t_023_FinalizeFilteredSnpsVCF.gcs_path
        File joint_snps_filtered_vcf_tbi = t_024_FinalizeFilteredSnpsTBI.gcs_path

        File joint_indels_filtered_vcf     = t_025_FinalizeFilteredIndelsVCF.gcs_path
        File joint_indels_filtered_vcf_tbi = t_026_FinalizeFilteredIndelsTBI.gcs_path

        File? annotated_joint_vcf     = annotated_vcf
        File? annotated_joint_vcf_tbi = annotated_vcf_tbi

        String? snpEff_summary = t_021_FinalizeSnpEffSummary.gcs_dir
        String? snpEff_genes = t_022_FinalizeSnpEffGenes.gcs_dir
    }
}

