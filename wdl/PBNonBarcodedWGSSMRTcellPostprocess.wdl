version 1.0

import "tasks/utils/AlignAndCheckFingerprintWGS.wdl" as major
import "tasks/Finalize.wdl" as FF
import "tasks/utils/GeneralUtils.wdl" as GU

workflow PBNonBarcodedWGSSMRTcellPostprocess {
    input {
        String movie_name
        String SM
        String LB

        File unaligned_bam
        File unaligned_pbi
        String experiment_type

        String fingerprint_store
        String sample_id_at_fpstore
        Boolean turn_off_fingperprint_check = false

        File ref_map_file

        String gcs_out_root_dir

        Boolean drop_per_base_N_pulse_tags = true
    }

    parameter_meta {
        fingerprint_store:   "Bucket name and prefix (gs://...) storing the fingerprinting VCFs"
        sample_id_at_fpstore: "A comma delimited string holding samples ID at fingerprint store"
        turn_off_fingperprint_check: "Please turn of fingerprint check if the reference is not GRCh38."

        ref_map_file:       "table indicating reference sequence and auxillary file locations"
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"

        experiment_type:    "type of experiment run (CLR, CCS, ISOSEQ, MASSEQ)"
        drop_per_base_N_pulse_tags: "when off, the puls tag (if available) in the input bam will be preserved; turning this one increases the output aligned bam."
    }

    Map[String, String] map_presets = {
        'CLR':    'SUBREAD',
        'CCS':    'CCS',
        'ISOSEQ': 'ISOSEQ',
        'MASSEQ': 'SUBREAD',
    }

    call major.AlignAndCheckFingerprintWGS {
        input:
            uBAM = unaligned_bam,
            uPBI = unaligned_pbi,
            bam_sample_name = SM,
            library = LB,
            experiment_type = experiment_type,

            turn_off_fingperprint_check = turn_off_fingperprint_check,
            fp_store = fingerprint_store,
            sample_id_at_store = sample_id_at_fpstore,
            ref_map_file = ref_map_file,
            drop_per_base_N_pulse_tags = drop_per_base_N_pulse_tags
    }

    ###################################################################################
    # Finalize data
    String workflow_name = "PBNonBarcodedWGSSMRTcellPostprocess"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name
    String outdir_aln     = outdir + '/alignments'
    String outdir_metrics = outdir + '/metrics'

    String mid =  if (experiment_type != "CLR") then "hifi" else "subreads"

    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = outdir_aln, file = AlignAndCheckFingerprintWGS.aligned_bam, name = movie_name + '.' + mid + '.bam' }
    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = outdir_aln, file = AlignAndCheckFingerprintWGS.aligned_bai, name = movie_name + '.' + mid + '.bai' }
    call FF.FinalizeToFile as FinalizeAlignedPbi { input: outdir = outdir_aln, file = AlignAndCheckFingerprintWGS.aligned_pbi, name = movie_name + '.' + mid + '.pbi' }

    call FF.FinalizeToFile as FinalizeAlnMetrics { input: outdir = outdir_metrics, file = AlignAndCheckFingerprintWGS.alignment_metrics_tar_gz }
    if (! turn_off_fingperprint_check) {
        call FF.FinalizeToFile as FinalizeFPDetails  { input: outdir = outdir_metrics, file = select_first([AlignAndCheckFingerprintWGS.fingerprint_detail_tar_gz]) }
    }

    ###################################################################################
    call GU.GetTodayDate as today {}

    output {
        # files for downstream WDLs
        File aligned_bam = FinalizeAlignedBam.gcs_path
        File aligned_bai = FinalizeAlignedBai.gcs_path
        File aligned_pbi = FinalizeAlignedPbi.gcs_path

        # metrics
        File aln_metric_tar_gz = FinalizeAlnMetrics.gcs_path

        String? fingerprint_check_result = AlignAndCheckFingerprintWGS.fp_status
        Float? fingerprint_check_LOD = AlignAndCheckFingerprintWGS.fp_lod_expected_sample
        File? fpcheck_tar_gz = FinalizeFPDetails.gcs_path

        ###############################
        String last_preprocessing_date = today.yyyy_mm_dd
    }
}
