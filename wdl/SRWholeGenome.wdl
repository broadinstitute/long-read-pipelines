version 1.0

######################################################################################
## A workflow that performs single sample variant calling on Illumina reads from
## one or more flow cells. The workflow merges multiple samples into a single BAM
## prior to variant calling.
######################################################################################

import "tasks/Utils.wdl" as Utils
import "tasks/SRUtils.wdl" as SRUTIL
import "tasks/CallVariantsIllumina.wdl" as VAR
import "tasks/HaplotypeCaller.wdl" as HC
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

    call COV.SampleLevelAlignedMetrics as coverage {
        input:
            aligned_bam = bam,
            aligned_bai = bai,
            ref_fasta   = ref_map['fasta'],
            bed_to_compute_coverage = bed_to_compute_coverage
    }

    String dir = outdir + "/alignments"

    call FF.FinalizeToFile as FinalizeBam { input: outdir = dir, file = bam, name = "~{participant_name}.bam" }
    call FF.FinalizeToFile as FinalizeBai { input: outdir = dir, file = bai, name = "~{participant_name}.bam.bai" }

    if (defined(bed_to_compute_coverage)) { call FF.FinalizeToFile as FinalizeRegionalCoverage { input: outdir = dir, file = select_first([coverage.bed_cov_summary]) } }

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

        call FF.FinalizeToFile as FinalizeDVPepperVcf { input: outdir = smalldir, file = select_first([CallVariantsWithDeepVariant.dvp_vcf])}
        call FF.FinalizeToFile as FinalizeDVPepperTbi { input: outdir = smalldir, file = select_first([CallVariantsWithDeepVariant.dvp_tbi])}
        call FF.FinalizeToFile as FinalizeDVPepperGVcf { input: outdir = smalldir, file = select_first([CallVariantsWithDeepVariant.dvp_g_vcf])}
        call FF.FinalizeToFile as FinalizeDVPepperGTbi { input: outdir = smalldir, file = select_first([CallVariantsWithDeepVariant.dvp_g_tbi])}
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

        call FF.FinalizeToFile as FinalizeHCVcf { input: outdir = smalldir, file = select_first([CallVariantsWithHaplotypeCaller.output_vcf])}
        call FF.FinalizeToFile as FinalizeHCTbi { input: outdir = smalldir, file = select_first([CallVariantsWithHaplotypeCaller.output_vcf_index])}
        call FF.FinalizeToFile as FinalizeHCGVcf { input: outdir = smalldir, file = select_first([CallVariantsWithHaplotypeCaller.output_gvcf])}
        call FF.FinalizeToFile as FinalizeHCGTbi { input: outdir = smalldir, file = select_first([CallVariantsWithHaplotypeCaller.output_gvcf_index])}
        call FF.FinalizeToFile as FinalizeHCBamOut { input: outdir = smalldir, file = select_first([CallVariantsWithHaplotypeCaller.bamout])}
        call FF.FinalizeToFile as FinalizeHCBaiOut { input: outdir = smalldir, file = select_first([CallVariantsWithHaplotypeCaller.bamout_index])}
    }

    output {
        File aligned_bam = FinalizeBam.gcs_path
        File aligned_bai = FinalizeBai.gcs_path

        Float aligned_num_reads = coverage.aligned_num_reads
        Float aligned_num_bases = coverage.aligned_num_bases
        Float aligned_frac_bases = coverage.aligned_frac_bases
        Float aligned_est_fold_cov = coverage.aligned_est_fold_cov

        Float aligned_read_length_mean = coverage.aligned_read_length_mean
        Float aligned_read_length_median = coverage.aligned_read_length_median
        Float aligned_read_length_stdev = coverage.aligned_read_length_stdev
        Float aligned_read_length_N50 = coverage.aligned_read_length_N50

        Float average_identity = coverage.average_identity
        Float median_identity = coverage.median_identity

        File? bed_cov_summary = FinalizeRegionalCoverage.gcs_path

        ########################################

        File? dvp_vcf   = FinalizeDVPepperVcf.gcs_path
        File? dvp_tbi   = FinalizeDVPepperTbi.gcs_path
        File? dvp_g_vcf = FinalizeDVPepperGVcf.gcs_path
        File? dvp_g_tbi = FinalizeDVPepperGTbi.gcs_path

        ########################################

        File? hc_vcf    = FinalizeHCVcf.gcs_path
        File? hc_tbi    = FinalizeHCTbi.gcs_path
        File? hc_g_vcf  = FinalizeHCGVcf.gcs_path
        File? hc_g_tbi  = FinalizeHCGTbi.gcs_path
        File? hc_bamout = FinalizeHCBamOut.gcs_path
        File? hc_baiout = FinalizeHCBaiOut.gcs_path
    }
}
