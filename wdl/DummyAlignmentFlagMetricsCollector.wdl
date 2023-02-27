version 1.0

import "tasks/metrics/AlignmentFlagMetrics.wdl" as worker

workflow DummyAlignmentFlagMetricsCollector {
    input {
        File aligned_bam
        File aligned_bai
    }

    call worker.AlignmentFlagMetrics as collector { input: aligned_bam = aligned_bam, aligned_bai = aligned_bai }
    output {
        Int aln_record_count = collector.aln_record_count
        Map [String, Float] flag_metrics = collector.flat_metrics
        Int summed_percentages = collector.summed_percentages
    }
}
