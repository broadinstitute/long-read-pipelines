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
import "tasks/HaplotypeCaller.wdl" as HC
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
        Int dvp_threads = 32
        Int dvp_memory = 128

        Int ploidy = 2

        Float snp_filter_level = 99.7
        Array[String] snp_recalibration_annotation_values = ["QD", "FS", "SOR", "MQRankSum", "ReadPosRankSum"]
        Array[Float] snp_recalibration_tranche_values = [100.0, 99.95, 99.9, 99.8, 99.6, 99.5, 99.4, 99.3, 99.0, 98.0, 97.0, 90.0 ]

        Float indel_filter_level = 99.0
        Array[String] indel_recalibration_annotation_values = ["QD", "FS", "SOR", "MQRankSum", "ReadPosRankSum"]
        Array[Float] indel_recalibration_tranche_values = [100.0, 99.95, 99.9, 99.5, 99.0, 97.0, 96.0, 95.0, 94.0, 93.5, 93.0, 92.0, 91.0, 90.0]

        File? bed_to_compute_coverage

        Array[String] contigs_names_to_ignore = ["RANDOM_PLACEHOLDER_VALUE"]  ## Required for ignoring any filtering - this is kind of a hack - TODO: fix the task!
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/SRWholeGenome/~{participant_name}"

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

    String bam_dir = outdir + "/alignments"

    call FF.FinalizeToFile as FinalizeBam { input: outdir = bam_dir, file = bam, name = "~{participant_name}.bam" }
    call FF.FinalizeToFile as FinalizeBai { input: outdir = bam_dir, file = bai, name = "~{participant_name}.bam.bai" }

    if (defined(bed_to_compute_coverage)) { call FF.FinalizeToFile as FinalizeRegionalCoverage { input: outdir = bam_dir, file = select_first([RegionalCoverage.cov_summary]) } }

    String metrics_dir = outdir + "/metrics"
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

    String smalldir = outdir + "/variants/small"

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
        call HC.CallVariantsWithHaplotypeCaller {
            input:
                bam               = bam,
                bai               = bai,
                sample_id         = participant_name,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                dbsnp_vcf         = ref_map["known_sites_vcf"],

                ploidy            = ploidy,

                prefix = participant_name + ".haplotype_caller",

                mito_contig = ref_map['mt_chr_name'],
                contigs_names_to_ignore = contigs_names_to_ignore,
        }

        call VARUTIL.IndelsVariantRecalibrator as TrainVQSROnHCIndelVariants {
            input:
                vcf = CallVariantsWithHaplotypeCaller.output_vcf,
                vcf_index = CallVariantsWithHaplotypeCaller.output_vcf_index,
                prefix = participant_name + ".indels",
                recalibration_tranche_values = snp_recalibration_tranche_values,
                recalibration_annotation_values = snp_recalibration_annotation_values,
                known_reference_variants = [ref_map["known_sites_vcf"]],
                known_reference_variants_index = [ref_map["known_sites_index"]],
                known_reference_variants_identifier = ["pfcrosses"],
                is_known = [true],
                is_training = [true],
                is_truth = [true],
                prior = [15],
                use_allele_specific_annotations = true,
                max_gaussians = 8,
        }

        call VARUTIL.SNPsVariantRecalibratorCreateModel as TrainVQSROnHCSnpVariants {
            input:
                vcf = CallVariantsWithHaplotypeCaller.output_vcf,
                vcf_index = CallVariantsWithHaplotypeCaller.output_vcf_index,
                prefix = participant_name + ".snps",
                recalibration_tranche_values = snp_recalibration_tranche_values,
                recalibration_annotation_values = snp_recalibration_annotation_values,
                known_reference_variants = [ref_map["known_sites_vcf"]],
                known_reference_variants_index = [ref_map["known_sites_index"]],
                known_reference_variants_identifier = ["pfcrosses"],
                is_known = [true],
                is_training = [true],
                is_truth = [true],
                prior = [15],
                use_allele_specific_annotations = true,
                max_gaussians = 8,
        }

        call VARUTIL.ApplyVqsr as ApplyVqsr {
            input:
                vcf = CallVariantsWithHaplotypeCaller.output_vcf,
                vcf_index = CallVariantsWithHaplotypeCaller.output_vcf_index,

                prefix = participant_name + ".vqsr_filtered",

                snps_recalibration = TrainVQSROnHCSnpVariants.recalibration,
                snps_recalibration_index = TrainVQSROnHCSnpVariants.recalibration_index,
                snps_tranches = TrainVQSROnHCSnpVariants.tranches,
                snp_filter_level = snp_filter_level,

                indels_recalibration = TrainVQSROnHCIndelVariants.recalibration,
                indels_recalibration_index = TrainVQSROnHCIndelVariants.recalibration_index,
                indels_tranches = TrainVQSROnHCIndelVariants.tranches,
                indel_filter_level = indel_filter_level,

                use_allele_specific_annotations = true,
        }

        call VARUTIL.RenameSingleSampleVcf as RenameSingleSampleVcf {
            input:
                vcf = ApplyVqsr.recalibrated_vcf,
                vcf_index = ApplyVqsr.recalibrated_vcf_index,
                prefix = participant_name + ".vqsr_filtered",
                new_sample_name = participant_name
        }

        # Create a Keyfile for finalization:
        File keyfile = ApplyVqsr.recalibrated_vcf_index

        # Finalize the raw Joint Calls:
        call FF.FinalizeToFile as FinalizeHCVcf { input: outdir = smalldir, keyfile = keyfile, file = CallVariantsWithHaplotypeCaller.output_vcf }
        call FF.FinalizeToFile as FinalizeHCTbi { input: outdir = smalldir, keyfile = keyfile, file = CallVariantsWithHaplotypeCaller.output_vcf_index }
        call FF.FinalizeToFile as FinalizeHCGVcf { input: outdir = smalldir, keyfile = keyfile, file = CallVariantsWithHaplotypeCaller.output_gvcf }
        call FF.FinalizeToFile as FinalizeHCGTbi { input: outdir = smalldir, keyfile = keyfile, file = CallVariantsWithHaplotypeCaller.output_gvcf_index }
        call FF.FinalizeToFile as FinalizeHCBamOut { input: outdir = smalldir, keyfile = keyfile, file = CallVariantsWithHaplotypeCaller.bamout }
        call FF.FinalizeToFile as FinalizeHCBaiOut { input: outdir = smalldir, keyfile = keyfile, file = CallVariantsWithHaplotypeCaller.bamout_index }

        # Finalize the VQSR files:
        call FF.FinalizeToFile as FinalizeIndelRecalFile { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCIndelVariants.recalibration }
        call FF.FinalizeToFile as FinalizeIndelRecalIndex { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCIndelVariants.recalibration_index }
        call FF.FinalizeToFile as FinalizeIndelRecalTranches { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCIndelVariants.tranches }
        call FF.FinalizeToFile as FinalizeIndelRecalModelReport { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCIndelVariants.model_report }

        call FF.FinalizeToFile as FinalizeSnpRecalFile { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.recalibration }
        call FF.FinalizeToFile as FinalizeSnpRecalIndex { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.recalibration_index }
        call FF.FinalizeToFile as FinalizeSnpRecalTranches { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.tranches }
        call FF.FinalizeToFile as FinalizeSnpRecalModelReport { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.model_report }

        # Finalize the reclibrated / filtered variants:
        call FF.FinalizeToFile as FinalizeHCVqsrVcf { input: outdir = smalldir, keyfile = keyfile, file = RenameSingleSampleVcf.new_sample_name_vcf }
        call FF.FinalizeToFile as FinalizeHCVqsrTbi { input: outdir = smalldir, keyfile = keyfile, file = RenameSingleSampleVcf.new_sample_name_vcf_index }
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

        File? bed_cov_summary = FinalizeRegionalCoverage.gcs_path

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
        File? hc_vqsr_vcf = FinalizeHCVqsrVcf.gcs_path
        File? hc_vqsr_tbi = FinalizeHCVqsrTbi.gcs_path

        File? vqsr_indel_recal_file         = FinalizeIndelRecalFile.gcs_path
        File? vqsr_indel_recal_file_index   = FinalizeIndelRecalIndex.gcs_path
        File? vqsr_indel_recal_tranches     = FinalizeIndelRecalTranches.gcs_path
        File? vqsr_indel_recal_model_report = FinalizeIndelRecalModelReport.gcs_path

        File? vqsr_snp_recal_file         = FinalizeSnpRecalFile.gcs_path
        File? vqsr_snp_recal_file_index   = FinalizeSnpRecalIndex.gcs_path
        File? vqsr_snp_recal_tranches     = FinalizeSnpRecalTranches.gcs_path
        File? vqsr_snp_recal_model_report = FinalizeSnpRecalModelReport.gcs_path
    }
}
