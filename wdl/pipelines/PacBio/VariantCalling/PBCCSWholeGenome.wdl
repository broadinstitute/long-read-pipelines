version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/VariantCalling/CallVariantsPBCCS.wdl" as VAR
import "../../../tasks/QC/SampleLevelAlignedMetrics.wdl" as COV

workflow PBCCSWholeGenome {

    meta {
        description: "A workflow that performs single sample variant calling on PacBio HiFi reads from one or more flow cells. The workflow merges multiple SMRT cells into a single BAM prior to variant calling."
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

        sites_vcf:     "for use with Clair"
        sites_vcf_tbi: "for use with Clair"

        run_dv_pepper_analysis:  "to turn on DV-Pepper analysis or not (non-trivial increase in cost and runtime)"
    }

    input {
        Array[File] aligned_bams

        File? bed_to_compute_coverage

        File ref_map_file

        String participant_name

        String gcs_out_root_dir

        Boolean call_svs = true
        Boolean fast_less_sensitive_sv = true

        Boolean call_small_variants = true
        Boolean call_small_vars_on_mitochondria = false
        File? sites_vcf
        File? sites_vcf_tbi

        Boolean run_dv_pepper_analysis = true
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCCSWholeGenome/" + participant_name
    String alignments_dir = outdir + "/alignments"

    # gather across (potential multiple) input CCS BAMs
    call Utils.MergeBams as MergeAllReads {
        input:
            bams = aligned_bams,
            outputBamName = "~{participant_name}.bam",
            outputBucket = alignments_dir,
            pacBioBams = true }

    call COV.SampleLevelAlignedMetrics as coverage {
        input:
            aligned_bam = MergeAllReads.merged_bam,
            aligned_bai = MergeAllReads.merged_bai
    }
    if (defined(bed_to_compute_coverage)) {
        call COV.MosDepthOverBed {
            input:
                bam = MergeAllReads.merged_bam,
                bai = MergeAllReads.merged_bai,
                bed = select_first([bed_to_compute_coverage]),
                output_bucket = alignments_dir
        }
    }

    ####################################################################################################
    if (call_svs || call_small_variants) {
        call VAR.CallVariants {
            input:
                bam               = MergeAllReads.merged_bam,
                bai               = MergeAllReads.merged_bai,
                prefix            = participant_name,
                output_bucket     = outdir,
                sample_id         = participant_name,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                regions_file      = ref_map[if (call_small_vars_on_mitochondria) then 'regions' else 'regionsNoM'],
                call_svs = call_svs,
                fast_less_sensitive_sv = fast_less_sensitive_sv,
                tandem_repeat_bed = ref_map['tandem_repeat_bed'],
                call_small_variants = call_small_variants,
                sites_vcf = sites_vcf,
                sites_vcf_tbi = sites_vcf_tbi,
                run_dv_pepper_analysis = run_dv_pepper_analysis
        }
    }

    output {
        File aligned_bam = MergeAllReads.merged_bam
        File aligned_bai = MergeAllReads.merged_bai
        File aligned_pbi = select_first([MergeAllReads.merged_pbi])

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
        File? pbsv_vcf = CallVariants.pbsv_vcf
        File? pbsv_tbi = CallVariants.pbsv_tbi

        File? sniffles_vcf = CallVariants.sniffles_vcf
        File? sniffles_tbi = CallVariants.sniffles_tbi

        File? clair_vcf = CallVariants.clair_vcf
        File? clair_tbi = CallVariants.clair_tbi

        File? clair_gvcf = CallVariants.clair_gvcf
        File? clair_gtbi = CallVariants.clair_gtbi

        File? dvp_vcf = CallVariants.dvp_vcf
        File? dvp_tbi = CallVariants.dvp_tbi
        File? dvp_g_vcf = CallVariants.dvp_g_vcf
        File? dvp_g_tbi = CallVariants.dvp_g_tbi
        File? dvp_phased_vcf = CallVariants.dvp_phased_vcf
        File? dvp_phased_tbi = CallVariants.dvp_phased_tbi
    }
}
