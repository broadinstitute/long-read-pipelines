version 1.0

import "tasks/Fusilli.wdl" as Fusilli

workflow FusilliCall {
    input {
        File ref_meta
        File ref_graph
        File ref_graph_colors

        Array[String] ref_ids
        Array[File] ref_paths

        Array[String] sample_ids
        Array[File] reads_fq1
        Array[File?] reads_fq2
        Array[File] cleaned_sample_graphs
        Array[File] cleaned_sample_graph_colors
        Array[Int] prune_thresholds

    }

    # Build links for each sample graph from both references and each sample's reads
    scatter(i in range(length(sample_ids))) {
        call Fusilli.ConstructSampleLinks as ConstructSampleLinks {
            input:
                sample_id = sample_ids[i],
                sample_graph = cleaned_sample_graphs[i],
                sample_graph_colors = cleaned_sample_graph_colors[i],

                ref_meta = ref_meta,
                ref_ids = ref_ids,
                ref_paths = ref_paths,

                reads_fq1 = reads_fq1[i],
                reads_fq2 = reads_fq2[i],

                prune_threshold = prune_thresholds[i]
        }
    }

    # The combined graph is a ccDBG combining all sample graphs, and is used to identify
    # variant k-mers in a later step below
    call Fusilli.BuildCombinedGraph as BuildCombinedGraph {
        input:
            sample_ids = sample_ids,
            cleaned_graphs = cleaned_sample_graphs,
            cleaned_graph_colors = cleaned_sample_graph_colors
    }

    call Fusilli.FindVariantKmers as FindVariantKmers {
        input:
            ref_graph = ref_graph,
            ref_graph_colors = ref_graph_colors,

            combined_graph = BuildCombinedGraph.combined_graph,
            combined_graph_colors = BuildCombinedGraph.combined_graph_colors
    }

    scatter(i in range(length(sample_ids))) {
        call Fusilli.AssembleVariantContigs as AssembleVariantContigs {
            input:
                sample_id = sample_ids[i],
                cleaned_graph = cleaned_sample_graphs[i],
                cleaned_graph_colors = cleaned_sample_graph_colors[i],
                linkdb = ConstructSampleLinks.sample_links[i],
                variant_kmers = FindVariantKmers.variant_kmers,
                nonref_kmers = FindVariantKmers.nonref_kmers
        }

    }
}
