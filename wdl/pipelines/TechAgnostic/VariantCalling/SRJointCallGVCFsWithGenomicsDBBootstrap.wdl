version 1.0

#############################################################################################################
## A workflow that performs joint calling on single-sample gVCFs from GATK4 HaplotypeCaller using GenomicsDB.
#############################################################################################################

import "../../../tasks/VariantCalling/SRJointGenotyping.wdl" as SRJOINT
import "../../../tasks/Utility/VariantUtils.wdl" as VARUTIL
import "../../../tasks/Utility/Utils.wdl" as UTILS
import "../../../tasks/Utility/Hail.wdl" as Hail
import "../../../tasks/TertiaryAnalysis/FunctionalAnnotation.wdl" as FUNK
import "../../../tasks/Utility/SGKit.wdl" as SGKit
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow SRJointCallGVCFsWithGenomicsDB {
    input {
        Array[File] gvcfs
        Array[File] gvcf_indices

        File ref_map_file

        File? interval_list

        Float heterozygosity = 0.001
        Float heterozygosity_stdev = 0.01
        Float indel_heterozygosity = 0.000125

        Array[File]?   annotation_bed_files
        Array[File]?   annotation_bed_file_indexes
        Array[String]? annotation_bed_file_annotation_names

        Array[String]? contig_list
        Array[File]? contig_interval_files

        File? snpeff_db
        String? genome_name 

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

    # Shard by contig for speed, defining contig_list/contig_interval_list to use
    Array[File] contig_interval_list_to_use = if defined (contig_interval_files) then select_first([contig_interval_files]) else t_001_MakeChrIntervalList.contig_interval_list_files
    
    scatter (idx_1 in range(length(contig_interval_list_to_use))) {

        String contig = if defined(contig_list) then select_first([contig_list])[idx_1] else t_001_MakeChrIntervalList.chrs[idx_1][0]
        File contig_interval_list = contig_interval_list_to_use[idx_1]

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
                runtime_attr_override = object {disk_gb: 10 + (3 * t_002_CreateSampleNameMap.total_gvcf_size_gb) + (2 * ceil(size(ref_map['fasta'], "GB"))), preemptible_tries: 0, max_retries:0}
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
                heterozygosity  = heterozygosity, 
                heterozygosity_stdev = heterozygosity_stdev,
                indel_heterozygosity = indel_heterozygosity,
                runtime_attr_override = object {preemptible_tries: 0} # Disable preemption for prototype.
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

    ##################################################################################
    #                                Hard Filtering                                  #
    ##################################################################################


    # Hard filter SNPs
    call VARUTIL.HardFilterVcfByGATKDefault_Snp as t_007_HardFilterSnps {
        input:
            vcf = t_006_MergeSitesOnlyVCFs.output_vcf,
            vcf_index = t_006_MergeSitesOnlyVCFs.output_vcf_index,
            prefix = prefix
    }

    # Hard filter indels
    call VARUTIL.HardFilterVcfByGATKDefault_Indel as t_008_HardFilterIndels {
        input:
            vcf = t_006_MergeSitesOnlyVCFs.output_vcf,
            vcf_index = t_006_MergeSitesOnlyVCFs.output_vcf_index,
            prefix = prefix
    }

    Array[File] filtered_vcfs = [t_007_HardFilterSnps.snp_filtered_vcf, t_008_HardFilterIndels.indel_filtered_vcf]
    Array[File] filtered_vcf_indices = [t_007_HardFilterSnps.snp_filtered_vcf_index, t_008_HardFilterIndels.indel_filtered_vcf_index]

    # Merge Snp and Indel VCF files
    call VARUTIL.MergeAndSortVCFsAllowOverlap as t_009_MergeSnpsAndIndels {
        input:
            vcfs = filtered_vcfs,
            vcf_indices = filtered_vcf_indices,
            ref_fasta_fai = ref_map["fai"],
            prefix = prefix
    }

    # Annotation by region
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

    # Functional annotation
    if (defined(snpeff_db)) {
        call FUNK.FunctionallyAnnotateVariants as t_011_FunctionallyAnnotate {
            input:
                    vcf = recalibrated_vcf,
                    snpeff_db = select_first([snpeff_db]),
                    genome_name = genome_name
        }
    }

    File finalized_vcf = select_first([t_011_FunctionallyAnnotate.annotated_vcf, recalibrated_vcf])
    File finalized_vcf_index = select_first([t_011_FunctionallyAnnotate.annotated_vcf_index, recalibrated_vcf_index])
    
    # Consolidate raw VCF files:
    call VARUTIL.GatherVcfs as t_012_GatherRawVcfs {
        input:
            input_vcfs = t_004_JointCallGVCFs.output_vcf,
            input_vcf_indices = t_004_JointCallGVCFs.output_vcf_index,
            prefix = prefix + ".raw.combined"
    }

    # Convert filtered SNPs to Zarr
    call SGKit.ConvertToZarrStore as t_013_ConvertToZarr {
        input:
            vcf = finalized_vcf,
            tbi = finalized_vcf_index,
            prefix = prefix
    }

    # Convert filtered SNPs to the output to a HAIL Matrix Table:
    call Hail.ConvertToHailMT as t_014_CreateHailMatrixTable {
        input:
            gvcf = finalized_vcf,
            tbi = finalized_vcf_index,
            reference = sub(sub(ref_map["fasta"], "^.*/", ""), "\.[fasta]*$", ""),
            ref_fasta = ref_map["fasta"],
            ref_fai = ref_map["fai"],
            prefix = prefix
    }

    ################################
    # Finalize the regular output files:
    ################################

    if (defined(gcs_out_root_dir)) {
        String concrete_gcs_out_root_dir = select_first([gcs_out_root_dir])
        String outdir = sub(concrete_gcs_out_root_dir, "/$", "") + "/SRJointCallGVCFsWithGenomicsDB/~{prefix}"

        String recalibration_dir = outdir + "/recalibration_files"
        String recalibration_model_dir = outdir + "/recalibration_files/model"
        String recalibration_results_dir = outdir + "/recalibration_files/results"
        String snpeff_results_dir = outdir + "/snpEff_results"

        File keyfile = t_014_CreateHailMatrixTable.completion_file

        call FF.FinalizeToDir as t_015_FinalizeGenomicsDB { input: outdir = outdir + "/GenomicsDB", keyfile = keyfile, files = t_003_ImportGVCFsIntoGenomicsDB.output_genomicsdb }

        call FF.FinalizeToFile as t_016_FinalizeRawVCF { input: outdir = outdir, keyfile = keyfile, file = t_012_GatherRawVcfs.output_vcf }
        call FF.FinalizeToFile as t_017_FinalizeRawTBI { input: outdir = outdir, keyfile = keyfile, file = t_012_GatherRawVcfs.output_vcf_index }

        call FF.FinalizeToFile as t_018_FinalizeFilteredVCF { input: outdir = outdir, keyfile = keyfile, file = t_009_MergeSnpsAndIndels.vcf }
        call FF.FinalizeToFile as t_019_FinalizeFilteredTBI { input: outdir = outdir, keyfile = keyfile, file = t_009_MergeSnpsAndIndels.tbi }

        if (defined(snpeff_db)) {
            call FF.FinalizeToFile as t_020_FinalizeSnpEffSummary { input: outdir = snpeff_results_dir, keyfile = keyfile, file = select_first([t_011_FunctionallyAnnotate.snpEff_summary]) }
            call FF.FinalizeToFile as t_021_FinalizeSnpEffGenes { input: outdir = snpeff_results_dir, keyfile = keyfile, file = select_first([t_011_FunctionallyAnnotate.snpEff_genes]) }
            
            call FF.FinalizeToFile as t_022_FinalizeAnnotatedFilteredVCF { input: outdir = outdir, keyfile = keyfile, file = select_first([t_011_FunctionallyAnnotate.annotated_vcf]) }
            call FF.FinalizeToFile as t_023_FinalizeAnnotatedFilteredTBI { input: outdir = outdir, keyfile = keyfile, file = select_first([t_011_FunctionallyAnnotate.annotated_vcf_index]) }
        }

        #####################################
        # Finalize the hard filtered files:
        #####################################
        call FF.FinalizeToFile as t_024_FinalizeFilteredSnpsVCF { input: outdir = outdir, keyfile = keyfile, file = t_007_HardFilterSnps.snp_filtered_vcf }
        call FF.FinalizeToFile as t_025_FinalizeFilteredSnpsTBI { input: outdir = outdir, keyfile = keyfile, file = t_007_HardFilterSnps.snp_filtered_vcf_index }

        call FF.FinalizeToFile as t_026_FinalizeFilteredIndelsVCF { input: outdir = outdir, keyfile = keyfile, file = t_008_HardFilterIndels.indel_filtered_vcf }
        call FF.FinalizeToFile as t_027_FinalizeFilteredIndelsTBI { input: outdir = outdir, keyfile = keyfile, file = t_008_HardFilterIndels.indel_filtered_vcf_index }

        #####################################
        # Finalize zarr/hail files:
        #####################################
        call FF.FinalizeToFile as t_028_FinalizeZarr { input: outdir = outdir, keyfile = keyfile, file = t_013_ConvertToZarr.zarr }
        call FF.FinalizeToFile as t_029_FinalizeHailMatrixTable {
            input:
                outdir = outdir,
                keyfile = keyfile,
                file = t_014_CreateHailMatrixTable.mt_tar
        }

        # Set up variable for outputs:
        Array[String] final_genomicsdb_location = [t_015_FinalizeGenomicsDB.gcs_dir]
    }

    # Make an alias for the functionally annotated data:
    if (defined(snpeff_db)) {
        File? annotated_vcf = if defined(gcs_out_root_dir) then select_first([t_022_FinalizeAnnotatedFilteredVCF.gcs_path]) else t_011_FunctionallyAnnotate.annotated_vcf
        File? annotated_vcf_tbi = if defined(gcs_out_root_dir) then select_first([t_023_FinalizeAnnotatedFilteredTBI.gcs_path]) else t_011_FunctionallyAnnotate.annotated_vcf_index

        File final_snpeff_summary = select_first([t_020_FinalizeSnpEffSummary.gcs_path, t_011_FunctionallyAnnotate.snpEff_summary])
        File final_snpeff_genes = select_first([t_021_FinalizeSnpEffGenes.gcs_path, t_011_FunctionallyAnnotate.snpEff_genes])
    }

    output {
        Array[String] genomicsDB = select_first([final_genomicsdb_location, t_003_ImportGVCFsIntoGenomicsDB.output_genomicsdb])

        File raw_joint_vcf     = select_first([t_016_FinalizeRawVCF.gcs_path, t_012_GatherRawVcfs.output_vcf])
        File raw_joint_vcf_tbi = select_first([t_017_FinalizeRawTBI.gcs_path, t_012_GatherRawVcfs.output_vcf_index])

        File joint_filtered_vcf     = select_first([t_018_FinalizeFilteredVCF.gcs_path, t_009_MergeSnpsAndIndels.vcf])
        File joint_filtered_vcf_tbi = select_first([t_019_FinalizeFilteredTBI.gcs_path, t_009_MergeSnpsAndIndels.tbi])

        File joint_zarr = select_first([t_028_FinalizeZarr.gcs_path, t_013_ConvertToZarr.zarr])
        File joint_mt = select_first([t_029_FinalizeHailMatrixTable.gcs_path, t_014_CreateHailMatrixTable.mt_tar])

        File joint_snps_filtered_vcf     = select_first([t_024_FinalizeFilteredSnpsVCF.gcs_path, t_007_HardFilterSnps.snp_filtered_vcf])
        File joint_snps_filtered_vcf_tbi = select_first([t_025_FinalizeFilteredSnpsTBI.gcs_path, t_007_HardFilterSnps.snp_filtered_vcf_index ])
        File joint_indels_filtered_vcf     = select_first([t_026_FinalizeFilteredIndelsVCF.gcs_path, t_008_HardFilterIndels.indel_filtered_vcf])
        File joint_indels_filtered_vcf_tbi = select_first([t_027_FinalizeFilteredIndelsTBI.gcs_path, t_008_HardFilterIndels.indel_filtered_vcf_index])

        File? annotated_joint_vcf     = annotated_vcf
        File? annotated_joint_vcf_tbi = annotated_vcf_tbi

        Array[File?] snpEff_summary = select_first([[final_snpeff_summary], []])
        Array[File?] snpEff_genes = select_first([[final_snpeff_genes], []])
    }
}

