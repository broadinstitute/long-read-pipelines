version 1.0

import "tasks/NanoPlot.wdl"

import "tasks/Utils.wdl"
import "tasks/utils/BAMutils.wdl"
import "tasks/utils/GeneralUtils.wdl" as GU

import "tasks/Finalize.wdl" as FF

workflow ProcessOnInstrumentDemuxedChunkRefFree {

    meta {
        desciption: "!!! WARN: THIS IS PROJECT-CENTER SPECIFIC !!! Given an on-instrument demultiplexed hifi_reads.bam, perform ref-independent prep work."
    }

    input {
        File uBAM

        String smrtcell_id
        String barcode_name

        String biosample_id
        String sample_id
        String downstream_sample_id

        String gcs_out_root_dir
    }

    ###################################################################################
    # prep work

    # where to store final results
    String workflow_name = "ProcessOnInstrumentDemuxedChunkRefFree"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name
    String outdir_ref_free = outdir + '/RefFree'

    ###################################################################################
    # stats
    call NanoPlot.NanoPlotFromUBam {input: bam = uBAM}

    call Utils.BamToFastq {input: bam = uBAM, prefix = "does_not_matter"}
    ###################################################################################
    # finalize
    String bc = barcode_name

    call BAMutils.GetReadGroupInfo as RG {input: uBAM = uBAM, keys = ['PU']}
    String movie_name = RG.read_group_info['PU']

    String bc_specific_fastq_out = outdir_ref_free + '/' + smrtcell_id + '/' + bc
    call FF.FinalizeToFile as FinalizeFQ { input: outdir = bc_specific_fastq_out, file = BamToFastq.reads_fq, name =  movie_name + '.' + bc + '.hifi.fq.gz' }

    ###################################################################################

    call GU.GetTodayDate as today {}

    output {

        File hifi_fq = FinalizeFQ.gcs_path
        Map[String, Float] hifi_stats = NanoPlotFromUBam.stats_map

        String last_postprocessing_date = today.yyyy_mm_dd
    }
}
