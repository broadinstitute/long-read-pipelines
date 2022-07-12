version 1.0

import "tasks/Fusilli.wdl" as Fusilli

workflow FusilliBuildAndCleanSampleGraph {
    input {
        String sample_id

        File reads_fq1
        File? reads_fq2
    }

    call Fusilli.BuildSampleGraph as BuildSampleGraph {
        input: sample_id = sample_id, reads_fq1 = reads_fq1, reads_fq2 = reads_fq2
    }

    Map[String, Float] clean_meta = read_json(BuildSampleGraph.clean_meta)

    output {
        File sample_graph = BuildSampleGraph.sample_graph
        File sample_graph_colors = BuildSampleGraph.sample_graph_colors
        File kmer_counts = BuildSampleGraph.kmer_counts
        File cleaned_graph = BuildSampleGraph.cleaned_graph
        File cleaned_graph_colors = BuildSampleGraph.cleaned_graph_colors
        Float est_coverage = clean_meta['est_coverage']
        Int prune_threshold = round(clean_meta['threshold'])
    }
}
