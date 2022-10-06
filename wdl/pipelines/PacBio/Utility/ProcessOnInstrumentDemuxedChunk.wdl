version 1.0

import "../../../tasks/Utility/PBUtils.wdl"
import "../../../tasks/Utility/BAMutils.wdl" as BU
import "../../../tasks/Utility/GeneralUtils.wdl" as GU

import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/Alignment/AlignAndCheckFingerprintCCS.wdl" as major
import "../../../tasks/QC/AlignedMetrics.wdl"

workflow ProcessOnInstrumentDemuxedChunk {

    meta {
        desciption: "!!! WARN: THIS IS PROJECT-CENTER SPECIFIC !!! Given an on-instrument demultiplexed hifi_reads.bam, perform alignment and QC check."
    }

    input {
        File uBAM

        String readgroup_id

        String bam_SM_field

        String fingerprint_store
        String sample_id_at_store
        Boolean turn_off_fingperprint_check = false

        File ref_map_file

        String gcs_out_root_dir
    }

    ###################################################################################
    # prep work

    # where to store final results
    String workflow_name = "ProcessOnInstrumentDemuxedChunk"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name

    ###################################################################################
    # generate PBI
    call PBUtils.PBIndex as Index {input: bam = uBAM}
    call BU.GetReadGroupInfo as RG {input: uBAM = uBAM, keys = ['SM', 'LB', 'PU']}

    # major work
    call major.AlignAndCheckFingerprintCCS {
        input:
            uBAM = uBAM,
            uPBI = Index.pbi,
            bam_sample_name = bam_SM_field,
            library = RG.read_group_info['LB'],

            turn_off_fingperprint_check = turn_off_fingperprint_check,
            fp_store = fingerprint_store,
            sample_id_at_store = sample_id_at_store,
            ref_map_file = ref_map_file
    }

    call AlignedMetrics.MosDepthWGS { input: bam = AlignAndCheckFingerprintCCS.aligned_bam, bai = AlignAndCheckFingerprintCCS.aligned_bai}
    ###################################################################################
    # finalize
    String movie_name = RG.read_group_info['PU']
    String bc_specific_aln_out    = outdir + '/alignments/' + readgroup_id
    String bc_specific_metric_out = outdir + "/metrics/"    + readgroup_id

    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = bc_specific_aln_out, file = AlignAndCheckFingerprintCCS.aligned_bam, name = readgroup_id + '.bam' }
    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = bc_specific_aln_out, file = AlignAndCheckFingerprintCCS.aligned_bai, name = readgroup_id + '.bai' }
    call FF.FinalizeToFile as FinalizeAlignedPbi { input: outdir = bc_specific_aln_out, file = AlignAndCheckFingerprintCCS.aligned_pbi, name = readgroup_id + '.pbi' }

    call FF.FinalizeToFile as FinalizeAlnMetrics { input: outdir = bc_specific_metric_out, file = AlignAndCheckFingerprintCCS.alignment_metrics_tar_gz }

    if (! turn_off_fingperprint_check) {
        call FF.FinalizeToFile as FinalizeFPDetails  { input: outdir = bc_specific_metric_out, file = select_first([AlignAndCheckFingerprintCCS.fingerprint_detail_tar_gz]) }
    }

    ###################################################################################

    call GU.GetTodayDate as today {}

    output {
        File aligned_bam = FinalizeAlignedBam.gcs_path
        File aligned_bai = FinalizeAlignedBai.gcs_path
        File aligned_pbi = FinalizeAlignedPbi.gcs_path
        Float wgs_cov = MosDepthWGS.wgs_cov

        Map[String, Float] alignment_metrics = AlignAndCheckFingerprintCCS.alignment_metrics
        File alignment_metrics_tar_gz = FinalizeAlnMetrics.gcs_path

        String movie = movie_name

        String? fingerprint_check_result = AlignAndCheckFingerprintCCS.fp_status
        Float? fingerprint_check_LOD = AlignAndCheckFingerprintCCS.fp_lod_expected_sample
        File? fingerprint_check_tar_gz = FinalizeFPDetails.gcs_path

        String last_processing_date = today.yyyy_mm_dd
    }
}
