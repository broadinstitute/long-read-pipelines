version 1.0

import "Structs.wdl"
import "NanoPlot.wdl" as NP
import "Finalize.wdl" as FF

workflow Figures {
    input {
        Array[File] summary_files

        String? gcs_output_dir
    }

    call NP.NanoPlotFromSummary as NanoPlotFromSummaryPDF { input: summary_files = summary_files, plot_type = "pdf" }
    call NP.NanoPlotFromSummary as NanoPlotFromSummaryPNG { input: summary_files = summary_files, plot_type = "png" }

    if (defined(gcs_output_dir)) {
        call FF.FinalizeToDir as FinalizeNanoPlotFromSummaryPDFs  { input: files = NanoPlotFromSummaryPDF.plots, outdir = gcs_output_dir + "/figures/pdf/" }
        call FF.FinalizeToDir as FinalizeNanoPlotFromSummaryPNGs  { input: files = NanoPlotFromSummaryPNG.plots, outdir = gcs_output_dir + "/figures/png/" }
        call FF.FinalizeToDir as FinalizeNanoPlotFromSummaryStats { input: files = [ NanoPlotFromSummaryPNG.stats ], outdir = gcs_output_dir + "/nanostats/" }
    }

    output {
        Array[File] NanoPlotFromSummaryPDFs = NanoPlotFromSummaryPDF.plots
        Array[File] NanoPlotFromSummaryPNGs = NanoPlotFromSummaryPNG.plots
        File NanoPlotFromSummaryStats = NanoPlotFromSummaryPNG.stats
    }
}

