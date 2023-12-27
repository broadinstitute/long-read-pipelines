version 1.0

import "../../../tasks/Utility/GeneralUtils.wdl" as GU

import "../../../tasks/Alignment/AlignHiFiUBAM.wdl" as ALN

import "../../TechAgnostic/Utility/AlignedBamQCandMetrics.wdl" as QCMetrics

import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ProcessOnInstrumentDemuxedChunk {

    meta {
        desciption: "Given an on-instrument demultiplexed hifi_reads.bam, perform alignment and QC check, and collect a few metrics."
    }

    parameter_meta {

        #########
        # inputs
        readgroup_id:
        "Unique ID for a readgroup, for example (D|E)A[0-9]{6}-<barcode>; no whitespaces allowed"

        bam_SM_field:
        "value to place in the SM field of the resulting BAM header's @RG line"

        platform:
        "PacBio platform used for generating the data; accepted value: [Sequel, Revio]"

        gcs_out_root_dir:
        "output files will be copied over there"

        qc_metrics_config_json:
        "A config json to for running the QC and metrics-collection sub-workflow 'AlignedBamQCandMetrics'"

        fingerprint_sample_id:
        "For fingerprint verification: the ID of the sample supposedly this BAM belongs to; note that the fingerprint VCF is assumed to be located at {fingerprint_store}/{fingerprint_sample_id}*.vcf(.gz)?"

        expected_sex_type:
        "If provided, triggers sex concordance check. Accepted value: [M, F, NA, na]"

        #########
        # outputs
        wgs_cov:
        "whole genome mean coverage"

        nanoplot_summ:
        "Summary on alignment metrics provided by Nanoplot (todo: study the value of this output)"

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
        String gcs_out_root_dir

        File  uBAM
        File? uPBI

        String readgroup_id
        String bam_SM_field

        String platform

        # args for optional QC subworkflows
        File? qc_metrics_config_json
        String? fingerprint_sample_id
        String? expected_sex_type

        File ref_map_file
        String disk_type = "SSD"
    }

    output {
        String last_processing_date = today.yyyy_mm_dd

        File aligned_bam = FinalizeAlignedBam.gcs_path
        File aligned_bai = FinalizeAlignedBai.gcs_path
        File aligned_pbi = FinalizeAlignedPbi.gcs_path

        String movie = AlignHiFiUBAM.movie

        # metrics block (caveat: always has to keep an eye on the QC subworkflow about outputs)
        Float wgs_cov                                   = QCandMetrics.wgs_cov
        Map[String, Float] nanoplot_summ                = QCandMetrics.nanoplot_summ
        Map[String, Float] sam_flag_stats               = QCandMetrics.sam_flag_stats

        # fingerprint
        Map[String, String]? fingerprint_check          = QCandMetrics.fingerprint_check

        # contam
        Float? contamination_est                        = QCandMetrics.contamination_est

        # sex concordance
        Map[String, String]? inferred_sex_info          = QCandMetrics.inferred_sex_info

        # methyl
        Map[String, String]? methyl_tag_simple_stats    = QCandMetrics.methyl_tag_simple_stats

        # file-based QC/metrics outputs all packed into a finalization map
        Map[String, String] aBAM_metrics_files          = QCandMetrics.aBAM_metrics_files
    }

    ###################################################################################
    # prep work

    # where to store final results
    String workflow_name = "ProcessOnInstrumentDemuxedChunk"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name

    # String bc_specific_out = outdir + '/' + readgroup_id
    String bc_specific_aln_out     = outdir + '/alignments/' + readgroup_id
    String bc_specific_metrics_out = outdir + "/metrics/"    + readgroup_id

    ###################################################################################
    # align
    call ALN.AlignHiFiUBAM as AlignHiFiUBAM { input:
        uBAM = uBAM,
        uPBI = uPBI,
        bam_sample_name = bam_SM_field,
        ref_map_file = ref_map_file,
        application = 'HIFI',
        disk_type = disk_type
    }

    File aBAM = AlignHiFiUBAM.aligned_bam
    File aBAI = AlignHiFiUBAM.aligned_bai
    File aPBI = AlignHiFiUBAM.aligned_pbi

    # save
    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = bc_specific_aln_out, file = aBAM, name = readgroup_id + '.bam' }
    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = bc_specific_aln_out, file = aBAI, name = readgroup_id + '.bai' }
    call FF.FinalizeToFile as FinalizeAlignedPbi { input: outdir = bc_specific_aln_out, file = aPBI, name = readgroup_id + '.pbi' }

    ###################################################################################
    # QC
    AlignedBamQCnMetricsConfig qcm_config = read_json(select_first([qc_metrics_config_json]))
    call QCMetrics.Work as QCandMetrics { input:
        bam = aBAM,
        bai = aBAI,

        tech = platform,

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

        output_prefix = readgroup_id, # String output_prefix = bam_sample_name + "." + flowcell
        gcs_out_root_dir = bc_specific_metrics_out,
    }

    ###################################################################################

    call GU.GetTodayDate as today {}
}
