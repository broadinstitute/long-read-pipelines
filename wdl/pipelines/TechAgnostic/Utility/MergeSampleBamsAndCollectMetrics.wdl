version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/PBUtils.wdl" as PB
import "../../../tasks/Utility/ONTUtils.wdl" as ONT
import "../../../tasks/Utility/BAMutils.wdl" as BU
import "../../../tasks/Utility/GeneralUtils.wdl" as GU

import "../../TechAgnostic/Utility/DystPeaker.wdl"    as QC2
import "../../TechAgnostic/Utility/FASTQstats.wdl"    as QC3

import "../../TechAgnostic/Utility/AlignedBamQCandMetrics.wdl" as QCMetrics
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow Work {

    meta {
        description: "Merge a sample's (potential) multiple SMRT-/flow-cells data, QC the BAM, and collect alignment metrics."
    }

    parameter_meta {

        #########
        # inputs
        sample_name:
        "value that should be in the SM field of the resulting BAM header's @RG line; this is used to validate that the input BAMs all have the expected value"

        tech:
        "LR technology used for generating the data; accepted value: [ONT, Sequel, Revio]"

        bams_suspected_to_contain_dup_record:
        "ONT-specific: sometimes, a single ONT flowcell's outputs have a (strange) read-duplication issue; if this issue wasn't explicitly fixed before the alignment step, we recommend turning this flag on here."

        qc_metrics_config_json:
        "A config json to for running the QC and metrics-collection sub-workflow 'AlignedBamQCandMetrics'"

        fingerprint_sample_id:
        "For fingerprint verification: the ID of the sample supposedly this BAM belongs to; note that the fingerprint VCF is assumed to be located at {fingerprint_store}/{fingerprint_sample_id}*.vcf(.gz)?"

        expected_sex_type:
        "If provided, triggers sex concordance check. Accepted value: [M, F, NA, na]"

        run_seqkit_stats:
        "if true, collect more metrics using seqkit; if processing HiFi data, consider getting this metrics via MergeHiFiFastQs; if alignning to multiple references, you just need to collect this once"

        short_reads_threshold:
        "if provided, will trigger read length metrics collection; a length threshold below which reads are classified as short; does not filter reads; if processing HiFi data, consider getting this metrics via MergeHiFiFastQs"

        gcs_out_root_dir:
        "output files will be copied over there"

        #########
        # outputs
        coverage:
        "whole genome mean coverage"

        # nanoplot_summ:
        # "Summary on alignment metrics provided by Nanoplot (todo: study the value of this output)"

        seqkit_stats:
        "A few metrics output by seqkit stats"

        read_len_summaries:
        "A few metrics summarizing the read length distribution"
        read_len_peaks:
        "Peaks of the read length distribution (heruistic)"
        read_len_deciles:
        "Deciles of the read length distribution"

        sam_flag_stats:
        "SAM flag stats"

        fingerprint_check:
        "Summary on (human) fingerprint checking results"

        contamination_est:
        "cross-(human)individual contamination estimation by VerifyBAMID2"

        inferred_sex_info:
        "Inferred sex concordance information if expected sex type is provided"

        methyl_tag_simple_stats:
        "Simple stats on the reads with & without SAM methylation tags (MM/ML)."

        aBAM_metrics_files:
        "A map where keys are summary-names and values are paths to files generated from the various QC/metrics tasks"
    }

    input {
        String gcs_out_dir

        # sample specific
        String sample_name
        Array[File] aligned_bams
        Array[File] aligned_bais

        String tech
        Boolean bams_suspected_to_contain_dup_record

        File? qc_metrics_config_json
        String? fingerprint_sample_id
        String? expected_sex_type

        Boolean run_seqkit_stats = false
        Int? short_reads_threshold

        File ref_map_file
        String disk_type = "SSD"
    }

    output {
        String last_processing_date = today.yyyy_mm_dd

        File aligned_bam = FinalizeBam.gcs_path
        File aligned_bai = FinalizeBai.gcs_path
        File? aligned_pbi = FinalizePbi.gcs_path

        # metrics block (caveat: always has to keep an eye on the QC subworkflow about outputs)
        Float coverage                               = QCandMetrics.wgs_cov
        # Map[String, Float] nanoplot_summ             = QCandMetrics.nanoplot_summ
        Map[String, Float] sam_flag_stats            = QCandMetrics.sam_flag_stats

        # metrics for ONT (for Hifi, this should be collected via the assembly route)
        Map[String, Float]? seqkit_stats        = FASTQstats.stats
        Map[String, String]? read_len_summaries = DystPeaker.read_len_summaries
        Array[Int]? read_len_peaks              = DystPeaker.read_len_peaks
        Array[Int]? read_len_deciles            = DystPeaker.read_len_deciles

        # fingerprint
        Map[String, String]? fingerprint_check       = QCandMetrics.fingerprint_check

        # contam
        Float? contamination_est                     = QCandMetrics.contamination_est

        # sex concordance
        Map[String, String]? inferred_sex_info       = QCandMetrics.inferred_sex_info

        # methyl
        Map[String, String]? methyl_tag_simple_stats = QCandMetrics.methyl_tag_simple_stats

        # file-based QC/metrics outputs all packed into a finalization map
        Map[String, String] aBAM_metrics_files       = metrics_files_out
    }

    ###################################################################################
    # prep work
    ###################################################################################

    ###########################################################
    # where to store final results
    # tree-strucuture
    # <gcs_out_dir>
    # ├── alignments/
    # │   ├── bam
    # │   ├── bai
    # │   ├── ...
    # │   ├── ...
    # │   ├── bam
    # │   ├── bai
    # │   └── metrics/
    # |   |   └── <sample_name>/
    # |   |   ...
    #         └── <sample_name>/
    String outdir = sub(gcs_out_dir, "/$", "") + "/alignments"
    String aln_metrics_out = sub(outdir, "/$", "") + "/metrics/~{sample_name}"

    ###########################################################
    # input validation
    scatter (pair in zip(aligned_bams, aligned_bais)) {
        call Utils.InferSampleName {input: bam = pair.left, bai = pair.right}
    }
    call Utils.CheckOnSamplenames {input: sample_names = InferSampleName.sample_name}
    if (InferSampleName.sample_name[0] != sample_name) {
        call Utils.StopWorkflow as SM_mismatch { input: reason = "Provided sample name and those encoded in the BAM(s) don't match."}
    }

    ###################################################################################
    # main responsibility: merge, dedup (ONT), index (PacBio)
    ###################################################################################
    # merge input BAMs, if necessary
    if (length(aligned_bams) > 1) {
        call BU.MergeBamsWithSamtools as MergeAllReads { input: bams = aligned_bams, out_prefix = sample_name }
    }

    File aBAM = select_first([MergeAllReads.merged_bam, aligned_bams[0]])
    File aBAI = select_first([MergeAllReads.merged_bai, aligned_bais[0]])

    ###########################################################
    # ont specific: sometimes there are literal duplicate reads
    if (('ONT'==tech) && bams_suspected_to_contain_dup_record) {
        call ONT.DeduplicateBam as RemoveONTDuplicates { input:
            aligned_bam = aBAM, aligned_bai = aBAI
        }
    }

    ###########################################################
    # save bam and index
    File sample_bam = select_first([RemoveONTDuplicates.corrected_bam, aBAM])
    File sample_bai = select_first([RemoveONTDuplicates.corrected_bai, aBAI])

    call FF.FinalizeToFile as FinalizeBam { input: outdir = outdir, file = sample_bam, name = "~{sample_name}.bam" }
    call FF.FinalizeToFile as FinalizeBai { input: outdir = outdir, file = sample_bai, name = "~{sample_name}.bam.bai" }

    ###########################################################
    # pacbio-specific index
    if ('ONT'!=tech) {
        call PB.PBIndex as PBIndexSampleReads { input: bam = sample_bam }
        call FF.FinalizeToFile as FinalizePbi { input: outdir = outdir, file = PBIndexSampleReads.pbi, name = "~{sample_name}.bam.pbi" }
    }

    ###################################################################################
    # QC/metrics
    ###################################################################################

    AlignedBamQCnMetricsConfig qcm_config = read_json(select_first([qc_metrics_config_json]))
    ###########################################################
    call QCMetrics.Work as QCandMetrics { input:
        bam = sample_bam,
        bai = sample_bai,

        tech = tech,

        cov_bed            = qcm_config.cov_bed,
        cov_bed_descriptor = qcm_config.cov_bed_descriptor,

        fingerprint_vcf_store = qcm_config.fingerprint_vcf_store,
        fingerprint_sample_id = fingerprint_sample_id,

        expected_sex_type = expected_sex_type,

        vbid2_config_json = qcm_config.vbid2_config_json,

        methyl_tag_check_bam_descriptor = qcm_config.methyl_tag_check_bam_descriptor,
        save_methyl_uncalled_reads = qcm_config.save_methyl_uncalled_reads,

        ref_map_file = ref_map_file,
        disk_type = disk_type,

        output_prefix = sample_name,
        gcs_out_root_dir = aln_metrics_out,
    }

    ##########
    if (defined(short_reads_threshold)) {
        call QC2.DystPeaker { input:
            input_file=sample_bam, input_is_bam=true,
            id=sample_name,
            short_reads_threshold=select_first([short_reads_threshold]),
            gcs_out_root_dir= aln_metrics_out
        }
        Map[String, String] read_len_out = {"read_len_hist": DystPeaker.read_len_hist,
                                            "raw_rl_file": DystPeaker.read_len_summaries['raw_rl_file']}
                                            # pack the outputs of standard, alignment-based metrics with more metrics
        call GU.MergeMaps as PackFilesOut { input:
            one = QCandMetrics.aBAM_metrics_files,
            two = read_len_out,
        }
    }

    ##########
    if (run_seqkit_stats) {
        call QC3.FASTQstats { input:
            reads = sample_bam,
            file_type = "BAM",
            seq_type = 'dna',  # we specifically turn off length filter here
        }
    }
    ##########
    Map[String, String] metrics_files_out = select_first([PackFilesOut.merged, QCandMetrics.aBAM_metrics_files])

    ###################################################################################
    call GU.GetTodayDate as today {}
}
