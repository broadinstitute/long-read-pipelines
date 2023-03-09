version 1.0

import "tasks/utils/AlignAndCheckFingerprintCCS.wdl" as major
import "tasks/utils/BAMutils.wdl"
import "tasks/Utils.wdl"
import "tasks/utils/GeneralUtils.wdl" as GU
import "tasks/Finalize.wdl" as FF

workflow ProcessRevioSMRTcell {
    meta {
        desciption:
        ""
    }

    parameter_meta {
        barcode_to_BAM_SM_ids:  ""
        barcode_to_ubams:       ""

        barcode_to_aligned_bam: ""
    }

    input {
        String smrtcell_id
        String movie_name

        Map[String, String] barcode_to_ubams

        Map[String, String] barcode_to_BAM_SM_ids

        String fingerprint_store
        Map[String, String] barcode_to_fingerprint_SM_ids
        Boolean turn_off_fingperprint_check = false

        File ref_map_file

        String gcs_out_root_dir
    }

    ###################################################################################
    # prep work

    # where to store final results
    String workflow_name = "ProcessRevioSMRTcell"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name
    String outdir_aln     = outdir + '/alignments'
    String outdir_metrics = outdir + '/metrics'

    ###################################################################################
    # major deal
    # alignment and alignment metrics collection
    Array[Pair[String, String]] unnecessary_coersion = barcode_to_ubams
    scatter (bc_n_bam in unnecessary_coersion) {
        String bc   = bc_n_bam.left
        String uBAM = bc_n_bam.right

        call BAMutils.GetReadGroupInfo as RG {input: uBAM = uBAM, keys = ['SM', 'LB']}

        call major.AlignAndCheckFingerprintCCS {
            input:
                uBAM = uBAM, # GetDemxedFilePaths.bam_path,
                uPBI = uBAM + ".pbi", # GetDemxedFilePaths.pbi_path,
                bam_sample_name = barcode_to_BAM_SM_ids[bc], # RG.read_group_info['SM'],

                turn_off_fingperprint_check = turn_off_fingperprint_check,
                fp_store = fingerprint_store,
                sample_id_at_store = barcode_to_fingerprint_SM_ids[bc],
                ref_map_file = ref_map_file
        }
        ###################################################################################
        # finalize each barcode
        String bc_specific_aln_out    = outdir + '/alignments/' + smrtcell_id + '/' + bc
        String bc_specific_metric_out = outdir + "/metrics/"    + smrtcell_id + '/' + bc

        call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = bc_specific_aln_out, file = AlignAndCheckFingerprintCCS.aligned_bam, name = movie_name + '.' + bc + '.bam' }
        call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = bc_specific_aln_out, file = AlignAndCheckFingerprintCCS.aligned_bai, name = movie_name + '.' + bc + '.bai' }
        # call FF.FinalizeToFile as FinalizeAlignedPbi { input: outdir = bc_specific_aln_out, file = AlignAndCheckFingerprintCCS.aligned_pbi, name = movie_name + '.' + bc + '.pbi' }

        call FF.FinalizeToFile as FinalizeAlnMetrics { input: outdir = bc_specific_metric_out, file = AlignAndCheckFingerprintCCS.alignment_metrics_tar_gz }

        # if (! turn_off_fingperprint_check) {
        #     call FF.FinalizeToFile as FinalizeFPDetails  { input: outdir = bc_specific_metric_out, file = select_first([AlignAndCheckFingerprintCCS.fingerprint_detail_tar_gz]) }
        # }
    }


    ###################################################################################
    # For Terra data tables
    call LocateBarcodeSpecificFoldersOrFiles as bc_2_aln_dir {input: barcodes_list = bc, finalized_dir_or_file_for_each_barcode = bc_specific_aln_out}
    # call FF.FinalizeToFile as finalize_bc_2_aln_dir_tsv {input: file = bc_2_aln_dir.barcode_2_gs_path, outdir = outdir_aln + '/' + smrtcell_id }

    call LocateBarcodeSpecificFoldersOrFiles as bc_2_aln_metric_dir {input: barcodes_list = bc, finalized_dir_or_file_for_each_barcode = FinalizeAlnMetrics.gcs_path}
    # call FF.FinalizeToFile as finalize_bc_2_aln_metrics_tsv {input: file = bc_2_aln_metric_dir.barcode_2_gs_path, outdir = outdir_metrics + '/' + smrtcell_id }
    # if (! turn_off_fingperprint_check) {
    #     call LocateBarcodeSpecificFoldersOrFiles as bc_2_fp_details_dir {input: barcodes_list = bc, finalized_dir_or_file_for_each_barcode = select_all(FinalizeFPDetails.gcs_path)}
    #     call FF.FinalizeToFile as finalize_bc_2_fp_details_tsv {input: file = bc_2_fp_details_dir.barcode_2_gs_path, outdir = outdir_metrics + '/' + smrtcell_id }
    # }

    # Array[String?] fp_res  = AlignAndCheckFingerprintCCS.fp_status
    # Array[Float?]  fp_lods = AlignAndCheckFingerprintCCS.fp_lod_expected_sample

    call GU.GetTodayDate as today {}

    output {
        Map[String, String] barcode_to_aln_files = read_map(bc_2_aln_dir.barcode_2_gs_path)

        # File barcode_2_aln_dir = finalize_bc_2_aln_dir_tsv.gcs_path

        # Array[String?] fingerprint_check_results = fp_res
        # Array[Float?] fingerprint_check_LODs = fp_lods

        Map[String, String] barcode_to_aln_metrics = read_map(bc_2_aln_metric_dir.barcode_2_gs_path)
        # File barcode_2_aln_metric_tar_gz = finalize_bc_2_aln_metrics_tsv.gcs_path
        # File? barcode_2_fpcheck_tar_gz = finalize_bc_2_fp_details_tsv.gcs_path

        String last_postprocessing_date = today.yyyy_mm_dd
    }
}

task LocateBarcodeSpecificFoldersOrFiles {
    meta {
        desciption: "Generates a 2-col TSV where the 1st col is the barcode name, and 2nd is the GCS \'folder\' or file for that barcode."
    }
    input {
        Array[String] barcodes_list
        Array[String] finalized_dir_or_file_for_each_barcode
    }

    command <<<
        set -eux
        paste <(cat ~{write_lines(barcodes_list)}) \
              <(cat ~{write_lines(finalized_dir_or_file_for_each_barcode)}) \
              > barcode_2_gs_path.tsv
    >>>

    output {
        File barcode_2_gs_path = "barcode_2_gs_path.tsv"
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
