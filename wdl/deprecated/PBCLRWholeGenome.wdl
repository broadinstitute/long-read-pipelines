version 1.0

import "../tasks/Utility/Utils.wdl" as Utils
import "../tasks/Utility/Finalize.wdl" as FF
import "../tasks/QC/SampleLevelAlignedMetrics.wdl" as COV

import "tasks/CallVariantsPBCLR.wdl" as VAR

workflow PBCLRWholeGenome {

    meta {
        description: "A workflow that performs single sample variant calling on PacBio CLR reads from one or more flow cells. The workflow merges multiple SMRT cells into a single BAM prior to variant calling."
    }
    parameter_meta {
        aligned_bams:       "GCS path to aligned BAM files"
        participant_name:   "name of the participant from whom these samples were obtained"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"

        call_svs:               "whether to call SVs"
        fast_less_sensitive_sv: "to trade less sensitive SV calling for faster speed"

        call_small_variants: "whether to call small variants"
        call_small_vars_on_mitochondria: "if false, will not attempt to call variants on mitochondria; if true, some samples might fail (caller feature) due to lack of signal"
    }

    input {
        Array[File] aligned_bams

        File? bed_to_compute_coverage

        File ref_map_file

        String participant_name

        String gcs_out_root_dir

        Boolean call_svs = true
        Boolean? fast_less_sensitive_sv = true

        Boolean call_small_variants = true
        Boolean? call_small_vars_on_mitochondria = false
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCLRWholeGenome/~{participant_name}"

    # gather across (potential multiple) input CLR BAMs
    String dir = outdir + "/alignments"
    call Utils.MergeBams as MergeAllReads {
        input:
            bams = aligned_bams,
            outputBamName = "~{participant_name}.bam",
            outputBucket = dir,
            checkSingleSample = true,
            pacBioBams = true
    }
    File bam = MergeAllReads.merged_bam
    File bai = MergeAllReads.merged_bai
    File pbi = select_first([MergeAllReads.merged_pbi])

    call COV.SampleLevelAlignedMetrics as coverage {
        input:
            aligned_bam = bam,
            aligned_bai = bai
    }
    if (defined(bed_to_compute_coverage)) {
        call COV.MosDepthOverBed {
            input:
                bam = bam,
                bai = bai,
                bed = select_first([bed_to_compute_coverage]),
                outputBucket = dir
        }
    }

    if (call_svs || call_small_variants) {

        # verify arguments are provided
        if (call_svs) {
            if (! defined(fast_less_sensitive_sv)) {call Utils.StopWorkflow as fast_less_sensitive_sv_not_provided {input: reason = "Calling SVs without specifying arg fast_less_sensitive_sv"}}
        }

        call VAR.CallVariants {
            input:
                bam               = bam,
                bai               = bai,

                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                tandem_repeat_bed = ref_map['tandem_repeat_bed'],

                prefix = participant_name,

                call_svs = call_svs,
                fast_less_sensitive_sv = select_first([fast_less_sensitive_sv]),

                call_small_variants = call_small_variants,
                call_small_vars_on_mitochondria = select_first([call_small_vars_on_mitochondria]),
        }

        String svdir = outdir + "/variants/sv"
        String smalldir = outdir + "/variants/small"

        if (call_svs) {
            call FF.FinalizeToFile as FinalizePBSV { input: outdir = svdir, file = select_first([CallVariants.pbsv_vcf]) }
            call FF.FinalizeToFile as FinalizePBSVtbi { input: outdir = svdir, file = select_first([CallVariants.pbsv_tbi]) }

            call FF.FinalizeToFile as FinalizeSniffles { input: outdir = svdir, file = select_first([CallVariants.sniffles_vcf]) }
            call FF.FinalizeToFile as FinalizeSnifflesTbi { input: outdir = svdir, file = select_first([CallVariants.sniffles_tbi]) }
        }
    }

    output {
        File aligned_bam = bam
        File aligned_bai = bai
        File aligned_pbi = pbi

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

        File? bed_cov_summary = MosDepthOverBed.regions
        ########################################
        File? pbsv_vcf = FinalizePBSV.gcs_path
        File? pbsv_tbi = FinalizePBSVtbi.gcs_path

        File? sniffles_vcf = FinalizeSniffles.gcs_path
        File? sniffles_tbi = FinalizeSnifflesTbi.gcs_path
    }
}
