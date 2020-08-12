version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/Structs.wdl"
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/NanoPlot.wdl" as NP
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/Finalize.wdl" as FF

workflow Figures {
    input {
        Array[File] summary_files

        String per
        String type
        String label

        String? gcs_output_dir
    }

    call NP.NanoPlotFromSummary as NanoPlotFromSummaryPDF { input: summary_files = summary_files, plot_type = "pdf" }
    call NP.NanoPlotFromSummary as NanoPlotFromSummaryPNG { input: summary_files = summary_files, plot_type = "png" }

    if (defined(gcs_output_dir)) {
        String plotdir  = sub(sub(gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://") + "/figures/~{per}_~{type}_~{label}"
        String statsdir = sub(sub(gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://") + "/metrics/~{per}_~{type}_~{label}"

        call FF.FinalizeToDir as FinalizeNanoPlotFromSummaryPDFs  { input: files = NanoPlotFromSummaryPDF.plots, outdir = plotdir + "/pdf/" }
        call FF.FinalizeToDir as FinalizeNanoPlotFromSummaryPNGs  { input: files = NanoPlotFromSummaryPNG.plots, outdir = plotdir + "/png/" }
        call FF.FinalizeToDir as FinalizeNanoPlotFromSummaryStats { input: files = [ NanoPlotFromSummaryPNG.stats ], outdir = statsdir + "/nanostats/" }
    }

    output {
        Array[File] NanoPlotFromSummaryPDFs = NanoPlotFromSummaryPDF.plots
        Array[File] NanoPlotFromSummaryPNGs = NanoPlotFromSummaryPNG.plots
        File NanoPlotFromSummaryStats = NanoPlotFromSummaryPNG.stats
    }
}

