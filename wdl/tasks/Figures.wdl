version 1.0

import "Structs.wdl"
import "NanoPlot.wdl" as NP
import "Finalize.wdl" as FF

workflow Figures {
    input {
        Array[File] summary_files

        String? gcs_output_dir
    }

    call NP.NanoPlotFromSummary { input: summary_files = summary_files }

    if (defined(gcs_output_dir)) {
        call FF.FinalizeToDir as FinalizeNanoPlotFromSummaryPNGs  { input: files = NanoPlotFromSummary.plots, outdir = gcs_output_dir + "/figures/png/" }
        call FF.FinalizeToDir as FinalizeNanoPlotFromSummaryStats { input: files = [ NanoPlotFromSummary.stats ], outdir = gcs_output_dir + "/nanostats/" }
    }

    output {
        Array[File] NanoPlotFromSummaryPNGs = NanoPlotFromSummary.plots
        File NanoPlotFromSummaryStats = NanoPlotFromSummary.stats
    }
}

