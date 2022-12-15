version 1.0

######################################################################################
## A workflow that performs single sample variant calling on Illumina reads from
## one or more flow cells. The workflow merges multiple samples into a single BAM
## prior to variant calling.
######################################################################################

import "tasks/SRUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/CallVariantsIllumina.wdl" as VAR
import "tasks/Finalize.wdl" as FF

import "tasks/SampleLevelAlignedMetrics.wdl" as COV

workflow SRWholeGenome {
    input {
        Array[File] aligned_bams
        Array[File] aligned_bais

        File? bed_to_compute_coverage

        File ref_map_file

        String participant_name

        String gcs_out_root_dir

        Boolean call_small_variants = true

        Boolean? run_dv_pepper_analysis = true
        Int? dvp_threads = 32
        Int? dvp_memory = 128
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
    if (call_small_variants) {

        # verify arguments are provided
        if (call_small_variants) {
            if (! defined(run_dv_pepper_analysis)) {call Utils.StopWorkflow as run_dv_pepper_analysis_not_provided {input: reason = "Unprovided arg run_dv_pepper_analysis"}}
            if (! defined(dvp_threads)) {call Utils.StopWorkflow as dvp_threads_not_provided {input: reason = "Unprovided arg dvp_threads"}}
        }

        call VAR.CallVariants {
            input:
                bam               = bam,
                bai               = bai,
                sample_id         = participant_name,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],

                prefix = participant_name,

                call_small_variants = call_small_variants,

                run_dv_pepper_analysis = select_first([run_dv_pepper_analysis]),
                dvp_threads = select_first([dvp_threads]),
                dvp_memory = select_first([dvp_memory]),
        }

        String smalldir = outdir + "/variants/small"

        if (call_small_variants) {
            if (select_first([run_dv_pepper_analysis])) {
                call FF.FinalizeToFile as FinalizeDVPepperVcf { input: outdir = smalldir, file = select_first([CallVariants.dvp_vcf])}
                call FF.FinalizeToFile as FinalizeDVPepperTbi { input: outdir = smalldir, file = select_first([CallVariants.dvp_tbi])}
                call FF.FinalizeToFile as FinalizeDVPepperGVcf { input: outdir = smalldir, file = select_first([CallVariants.dvp_g_vcf])}
                call FF.FinalizeToFile as FinalizeDVPepperGTbi { input: outdir = smalldir, file = select_first([CallVariants.dvp_g_tbi])}
            }
        }
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

        File? dvp_vcf = FinalizeDVPepperVcf.gcs_path
        File? dvp_tbi = FinalizeDVPepperTbi.gcs_path
        File? dvp_g_vcf = FinalizeDVPepperGVcf.gcs_path
        File? dvp_g_tbi = FinalizeDVPepperGTbi.gcs_path
    }
}
