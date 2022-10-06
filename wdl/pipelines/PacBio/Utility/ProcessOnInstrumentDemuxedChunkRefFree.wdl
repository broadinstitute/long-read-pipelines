version 1.0

import "../../../tasks/Visualization/NanoPlot.wdl"

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/BAMutils.wdl" as BU
import "../../../tasks/Utility/GeneralUtils.wdl" as GU

import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ProcessOnInstrumentDemuxedChunkRefFree {

    meta {
        desciption: "!!! WARN: THIS IS PROJECT-CENTER SPECIFIC !!! Given an on-instrument demultiplexed hifi_reads.bam, perform ref-independent prep work."
    }

    input {
        File uBAM

        String readgroup_id

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
    call NanoPlot.NanoPlotFromUBam {input: uBAM = uBAM}

    call Utils.BamToFastq {input: bam = uBAM, prefix = "does_not_matter"}
    ###################################################################################
    # finalize
    call BU.GetReadGroupInfo as RG {input: uBAM = uBAM, keys = ['PU']}
    String movie_name = RG.read_group_info['PU']

    String bc_specific_fastq_out = outdir_ref_free + '/' + readgroup_id
    call FF.FinalizeToFile as FinalizeFQ { input: outdir = bc_specific_fastq_out, file = BamToFastq.reads_fq, name =  readgroup_id + '.hifi.fq.gz' }

    ###################################################################################

    call GU.GetTodayDate as today {}

    output {
        String movie = movie_name

        File hifi_fq = FinalizeFQ.gcs_path
        Map[String, Float] hifi_stats = NanoPlotFromUBam.stats_map

        String last_processing_date = today.yyyy_mm_dd
    }
}
