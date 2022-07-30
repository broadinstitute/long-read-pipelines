version 1.0

import "tasks/NanoPlot.wdl"

workflow HifiUstats {

    input {
        File unaligned_file
        Boolean file_is_bam

    }

    String format = if file_is_bam then "ubam" else "fastq"
    call NanoPlot.NanoPlotFromUnAligned as Metrics {input: unaligned_file = unaligned_file, format=format}

    output {
        Map[String, Float] hifi_stats_map = Metrics.stats_map
    }
}