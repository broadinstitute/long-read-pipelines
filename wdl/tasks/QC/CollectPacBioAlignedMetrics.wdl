version 1.0

import "../Utility/PBUtils.wdl" as PB
import "../Visualization/NanoPlot.wdl" as NP

workflow CollectPacBioAlignedMetrics {

    meta {
        desciption:
        "Collect a few custom metrics on the alignments."
    }
    parameter_meta {
        custom_aln_metrics_summary: "A 2-col TSV holding custom metrics on the alignments: 1-col is attribute name, 2-col is value."
        nanoplot_stats: "Nanoplot stats file (NanoStats.txt) on the input BAM."
        nanoplot_pngs: "Nanoplot figures on the input BAM."
    }

    input {
        File aligned_bam
        File aligned_bai
        File aligned_pbi
    }

    # note: this may not matter anymore if the input bam is HiFi only?
    call PB.SummarizePBI as SummarizeAlignedPBI    { input: pbi = aligned_pbi }
    call PB.SummarizePBI as SummarizeAlignedQ5PBI  { input: pbi = aligned_pbi, qual_threshold = 5 }
    call PB.SummarizePBI as SummarizeAlignedQ7PBI  { input: pbi = aligned_pbi, qual_threshold = 7 }
    call PB.SummarizePBI as SummarizeAlignedQ10PBI { input: pbi = aligned_pbi, qual_threshold = 10 }
    call PB.SummarizePBI as SummarizeAlignedQ12PBI { input: pbi = aligned_pbi, qual_threshold = 12 }
    call PB.SummarizePBI as SummarizeAlignedQ15PBI { input: pbi = aligned_pbi, qual_threshold = 15 }

    call NP.NanoPlotFromBam { input: bam = aligned_bam, bai = aligned_bai }

    call CustomMetricsSummaryToFile {
        input:
        attributes = ["num_reads_Q5", "num_reads_Q7", "num_reads_Q10", "num_reads_Q12", "num_reads_Q15",
                      "aligned_num_reads", "aligned_num_bases", "aligned_frac_bases",
                      "aligned_read_length_mean", "aligned_read_length_median", "aligned_read_length_stdev", "aligned_read_length_N50",
                      "average_identity", "median_identity"],
        values = [SummarizeAlignedQ5PBI.results['reads'], SummarizeAlignedQ7PBI.results['reads'], SummarizeAlignedQ10PBI.results['reads'], SummarizeAlignedQ12PBI.results['reads'], SummarizeAlignedQ15PBI.results['reads'],
                  NanoPlotFromBam.stats_map['number_of_reads'], NanoPlotFromBam.stats_map['number_of_bases_aligned'], NanoPlotFromBam.stats_map['fraction_bases_aligned'],
                  NanoPlotFromBam.stats_map['mean_read_length'], NanoPlotFromBam.stats_map['median_read_length'], NanoPlotFromBam.stats_map['read_length_stdev'], NanoPlotFromBam.stats_map['n50'],
                  NanoPlotFromBam.stats_map['average_identity'], NanoPlotFromBam.stats_map['median_identity']]
    }

    output {
        File custom_aln_metrics_summary = CustomMetricsSummaryToFile.custom_aln_metrics_summary
        File nanoplot_stats = NanoPlotFromBam.stats
        Map[String, Float] alignment_metrics = NanoPlotFromBam.stats_map
        Array[File] nanoplot_pngs = NanoPlotFromBam.plots
    }
}

task CustomMetricsSummaryToFile {
    meta {
        desciption:
        "Format a few custom metrics on the alignments into a 2-col TSV: 1-col is attribute name, 2-col is value."
    }

    input {
        Array[String] attributes
        Array[String] values
    }

    command <<<
        set -eux
        paste ~{write_lines(attributes)} ~{write_lines(values)} > "alignment.metrics.tsv"
    >>>

    output {
        File custom_aln_metrics_summary = "alignment.metrics.tsv"
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
