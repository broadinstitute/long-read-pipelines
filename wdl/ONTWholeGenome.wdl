version 1.0

######################################################################################
## A workflow that performs single sample variant calling on Oxford Nanopore reads
## from one or more flow cells. The workflow merges multiple samples into a single BAM
## prior to variant calling.
######################################################################################

import "tasks/ONTUtils.wdl" as ONT
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF

import "tasks/SampleLevelAlignedMetrics.wdl" as COV

workflow ONTWholeGenome {
    input {
        Array[File] aligned_bams
        Array[File] aligned_bais
        Boolean bams_suspected_to_contain_dup_record

        File? bed_to_compute_coverage

        File ref_map_file

        String participant_name

        String gcs_out_root_dir

        Boolean call_svs = true
        Boolean? fast_less_sensitive_sv = true

        Boolean call_small_variants = true
        Boolean? call_small_vars_on_mitochondria = false
        File? sites_vcf
        File? sites_vcf_tbi

        Boolean? run_dv_pepper_analysis = true
        Int? dvp_threads = 32
        Int? dvp_memory = 128
        File? ref_scatter_interval_list_locator
        File? ref_scatter_interval_list_ids
    }

    parameter_meta {
        aligned_bams:       "GCS path to aligned BAM files"
        aligned_bais:       "GCS path to aligned BAM file indices"
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
        ref_scatter_interval_list_locator: "A file holding paths to interval_list files; needed only when running DV-Pepper"
        ref_scatter_interval_list_ids:     "A file that gives short IDs to the interval_list files; needed only when running DV-Pepper"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTWholeGenome/~{participant_name}"

    # gather across (potential multiple) input raw BAMs
    if (length(aligned_bams) > 1) {
        scatter (pair in zip(aligned_bams, aligned_bais)) {
            call Utils.InferSampleName {input: bam = pair.left, bai = pair.right}
        }
        call Utils.CheckOnSamplenames {input: sample_names = InferSampleName.sample_name}

        call Utils.MergeBams as MergeAllReads { input: bams = aligned_bams, prefix = participant_name }
    }

    File bam = select_first([MergeAllReads.merged_bam, aligned_bams[0]])
    File bai = select_first([MergeAllReads.merged_bai, aligned_bais[0]])

    if (bams_suspected_to_contain_dup_record) {
        call Utils.DeduplicateBam as RemoveDuplicates {
            input: aligned_bam = bam, aligned_bai = bai
        }
    }
    File usable_bam = select_first([RemoveDuplicates.corrected_bam, bam])
    File usable_bai = select_first([RemoveDuplicates.corrected_bai, bai])

    #call COV.SampleLevelAlignedMetrics as coverage {
    #    input:
    #        aligned_bam = usable_bam,
    #        aligned_bai = usable_bai,
    #        ref_fasta   = ref_map['fasta'],
    #        bed_to_compute_coverage = bed_to_compute_coverage
    #}

    String dir = outdir + "/alignments"

    call FF.FinalizeToFile as FinalizeBam { input: outdir = dir, file = usable_bam, name = "~{participant_name}.bam" }
    call FF.FinalizeToFile as FinalizeBai { input: outdir = dir, file = usable_bai, name = "~{participant_name}.bam.bai" }

    output {
        File merged_bam = FinalizeBam.gcs_path
        File merged_bai = FinalizeBai.gcs_path

        #Float aligned_num_reads = coverage.aligned_num_reads
        #Float aligned_num_bases = coverage.aligned_num_bases
        #Float aligned_frac_bases = coverage.aligned_frac_bases
        #Float aligned_est_fold_cov = coverage.aligned_est_fold_cov

        #Float aligned_read_length_mean = coverage.aligned_read_length_mean
        #Float aligned_read_length_median = coverage.aligned_read_length_median
        #Float aligned_read_length_stdev = coverage.aligned_read_length_stdev
        #Float aligned_read_length_N50 = coverage.aligned_read_length_N50

        #Float average_identity = coverage.average_identity
        #Float median_identity = coverage.median_identity

        ########################################
    }
}
