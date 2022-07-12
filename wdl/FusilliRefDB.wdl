version 1.0

import "tasks/Fusilli.wdl" as Fusilli

workflow FusilliRefDB {
    input {
        Array[String] ref_ids
        Array[File] ref_paths

        String gcs_output_dir
    }

    call Fusilli.BuildGraph as BuildGraph {
        input: ref_ids = ref_ids, ref_paths = ref_paths
    }

    scatter(ref in zip(ref_ids, ref_paths)) {
        call Fusilli.BuildRefLinks as BuildRefLinks {
            input: graph = BuildGraph.graph, graph_colors = BuildGraph.graph_colors, ref_meta = BuildGraph.ref_meta,
                   ref_id = ref.left, ref_path = ref.right
        }
    }

    call Fusilli.FinalizeDB as Finalize {
        input:
            graph = BuildGraph.graph,
            graph_colors = BuildGraph.graph_colors,
            ref_meta = BuildGraph.ref_meta,

            ref_ids = ref_ids,
            ref_links = BuildRefLinks.ref_links,
            mm2_indices = BuildRefLinks.minimap2_index,

            gcs_output_dir = gcs_output_dir
    }

    output {
        String ref_graph = Finalize.ref_graph
        String ref_graph_colors = Finalize.ref_graph_colors
    }
}



