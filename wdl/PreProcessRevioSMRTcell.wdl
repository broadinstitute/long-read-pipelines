version 1.0

import "tasks/CCSLima.wdl"
import "tasks/SMRTtools.wdl"
import "tasks/metrics/CollectSMRTCellUnalignedMetrics.wdl" as uBAMCustomMetrics
import "tasks/PBUtils.wdl" as PB
import "tasks/Finalize.wdl" as FF
import "tasks/Utils.wdl"
import "tasks/utils/GeneralUtils.wdl" as GU


import "tasks/ParseRevioFolder.wdl" as Parser


workflow PreProcessRevioSMRTcell {
    meta {
        desciption:
        "A workflow for preprocessing barcoded (potentially multiplexed) Revio SMRT cell. The whole data folder is assumed to be mirrored onto a cloud bucket."
    }
    input {
        String smrtcell_folder

        String gcs_out_root_dir
    }
    parameter_meta {
        smrtcell_folder: "The (cloud) location of SMRT cell folder"
    }

    String workflow_name = "PreProcessRevioSMRTcell"

    call Parser.ParseRevioFolder { input: smrtcell_folder = smrtcell_folder}

    # todo
    # parse the on-instrument metrics and pick to output
    # run qc report from smrttools
    # run barcode report from smrttools
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name + "/" + ParseRevioFolder.movie
    # String outdir_metrics = outdir + "/metrics"

    Array[Pair[String, String]] coerced_pair_from_map = ParseRevioFolder.bc_2_bam
    scatter (pair in coerced_pair_from_map) {
        String barcode = pair.left
        String hifi_ubam = pair.right
        call Utils.BamToFastq { input: bam = hifi_ubam, prefix=basename(hifi_ubam, ".bam")}
        call FF.FinalizeToFile as FinalizeFq { input:  outdir = outdir, file = BamToFastq.reads_fq, }
    }
    call Parser.CoercePairsOfArrayToMap as AssociateToFastQ {
        input: keys = barcode, values = FinalizeFq.gcs_path
    }

    call GU.GetTodayDate as today {}

    output {
        String last_preprocessing_date = today.yyyy_mm_dd
        String run_id = ParseRevioFolder.run_id
        String cell_index = ParseRevioFolder.cell_idex

        String movie = ParseRevioFolder.movie

        Map[String, String] barcode_2_bam = ParseRevioFolder.bc_2_bam
        Map[String, String] barcode_2_fq  = AssociateToFastQ.res

        Map[String, String] barcode_2_aliquot_id   = ParseRevioFolder.bc_2_aliquot_id
        Map[String, String] barcode_2_biosample_id = ParseRevioFolder.bc_2_biosample_id
        Map[String, String] barcode_2_lims_id      = ParseRevioFolder.bc_2_lims_id
    }
}


