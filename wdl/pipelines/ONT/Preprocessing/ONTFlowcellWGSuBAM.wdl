version 1.0

import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/Alignment/AlignONTWGSuBAM.wdl" as Alignment

import "../../TechAgnostic/Utility/AlignedBamQCandMetrics.wdl" as QCMetrics

import "../../TechAgnostic/Utility/DystPeaker.wdl"    as QC2
import "../../TechAgnostic/Utility/FASTQstats.wdl"    as QC3

workflow ONTFlowcellWGSuBAM {

    meta {
        description:
        "Align one ONT flowcell's reads to a reference genome. Note this is NOT suitable for flowcells that have multiple basecall directories."
    }

    parameter_meta {
        # input
        uBAM:      "GCS path to the unaligned BAM of the flowcell. Note that it's assumed the infamous duplicated read issue has been handled."
        flowcell:  "flowcell ID"
        bam_SM_ID: "sample name to use in the aligned BAM output."
        uBAM_tags_to_preserve: "SAM tags to transfer to the alignment bam, from the input BAM. Please be careful with this: methylation tags MM, ML are usually what you want to keep at minimum."

        longer_ont_read_hint:
        "hint that the input reads are longer (~>20kb N50) ONT reads; this is used to bump the memory to the alignment task"

        run_seqkit_stats:
        "if true, collect more metrics using seqkit; if alignning to multiple references, you just need to collect this once"

        short_reads_threshold:
        "if provided, will trigger read length metrics collection; a length threshold below which reads are classified as short; does not filter reads"

        aln_disk_type: "An optimization specifying which type of disk to use for the minimap2 alignment step."

        ref_map_file:       "table indicating reference sequence and auxillary file locations"
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"

        #########
        # outputs
        wgs_cov:
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
        File uBAM
        String flowcell
        String bam_SM_ID
        Boolean longer_read_hint = false
        Array[String] uBAM_tags_to_preserve = ['MM', 'ML']

        # args for optional QC subworkflows
        File? qc_metrics_config_json
        String? fingerprint_sample_id
        String? expected_sex_type

        Boolean run_seqkit_stats = false
        Int? short_reads_threshold

        File ref_map_file
        String gcs_out_root_dir
        String aln_disk_type = 'SSD'
    }

    output {
        String last_processing_date = today.yyyy_mm_dd

        File aligned_bam = FinalizeAlignedBam.gcs_path
        File aligned_bai = FinalizeAlignedBai.gcs_path

        Float wgs_cov                                   = QCandMetrics.wgs_cov
        # Map[String, Float] nanoplot_summ                = QCandMetrics.nanoplot_summ
        Map[String, Float] sam_flag_stats               = QCandMetrics.sam_flag_stats

        Map[String, Float]? seqkit_stats = FASTQstats.stats

        Map[String, String]? read_len_summaries = DystPeaker.read_len_summaries
        Array[Int]? read_len_peaks = DystPeaker.read_len_peaks
        Array[Int]? read_len_deciles = DystPeaker.read_len_deciles

        # fingerprint
        Map[String, String]? fingerprint_check          = QCandMetrics.fingerprint_check

        # contam
        Float? contamination_est                        = QCandMetrics.contamination_est

        # sex concordance
        Map[String, String]? inferred_sex_info          = QCandMetrics.inferred_sex_info

        # methyl
        Map[String, String]? methyl_tag_simple_stats    = QCandMetrics.methyl_tag_simple_stats

        # file-based QC/metrics outputs all packed into a finalization map
        Map[String, String] aBAM_metrics_files          = files_out
    }

    ###################################################################################
    # prep work
    String workflow_name = "ONTFlowcellWGSuBAM"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_name}/~{flowcell}"
    String outdir_aln     = outdir + '/alignments'
    String outdir_metrics = outdir + '/metrics'

    ###################################################################################
    # align & save
    call Alignment.AlignONTWGSuBAM as ALN {
        input:
            uBAM = uBAM,
            longer_read_hint = longer_read_hint,
            uBAM_tags_to_preserve = uBAM_tags_to_preserve,
            bam_sample_name = bam_SM_ID,
            flowcell = flowcell,
            ref_map_file = ref_map_file,
            aln_disk_type = aln_disk_type
    }
    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = outdir_aln, file = ALN.aligned_bam }
    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = outdir_aln, file = ALN.aligned_bai }

    ###################################################################################
    # QC
    String readgroup_id = bam_SM_ID + "." + flowcell
    AlignedBamQCnMetricsConfig qcm_config = read_json(select_first([qc_metrics_config_json]))
    call QCMetrics.Work as QCandMetrics { input:
        bam = ALN.aligned_bam,
        bai = ALN.aligned_bai,

        tech = "ONT",

        cov_bed            = qcm_config.cov_bed,
        cov_bed_descriptor = qcm_config.cov_bed_descriptor,

        fingerprint_vcf_store = qcm_config.fingerprint_vcf_store,
        fingerprint_sample_id = fingerprint_sample_id,

        expected_sex_type = expected_sex_type,

        vbid2_config_json = qcm_config.vbid2_config_json,

        methyl_tag_check_bam_descriptor = qcm_config.methyl_tag_check_bam_descriptor,
        save_methyl_uncalled_reads = qcm_config.save_methyl_uncalled_reads,

        ref_map_file = ref_map_file,
        disk_type = aln_disk_type,

        output_prefix = readgroup_id,
        gcs_out_root_dir = outdir_metrics,
    }

    ##########
    if (defined(short_reads_threshold)) {
        call QC2.DystPeaker { input:
            input_file=uBAM, input_is_bam=true,
            id=readgroup_id,
            short_reads_threshold=select_first([short_reads_threshold]),
            gcs_out_root_dir= outdir_metrics
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
            reads = uBAM,
            file_type = "BAM",
            seq_type = 'dna',  # we specifically turn off length filter here
        }
    }

    ##########
    Map[String, String] files_out = select_first([PackFilesOut.merged, QCandMetrics.aBAM_metrics_files])

    ###################################################################################
    call GU.GetTodayDate as today {}
}
