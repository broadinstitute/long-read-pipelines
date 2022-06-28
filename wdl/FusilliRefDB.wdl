version 1.0

import "tasks/Fusilli.wdl" as Fusilli
import "tasks/Finalize.wdl" as FF

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
        call Fusilli.BuildLinks as BuildLinks {
            input: graph = BuildGraph.graph, graph_colors = BuildGraph.graph_colors, ref_meta = BuildGraph.ref_meta,
                   ref_id = ref[0], ref_path = ref[1]
        }
    }

    call Fusilli.FinalizeDB as Finalize {
        input:
            graph = BuildGraph.graph,
            graph_colors = BuildGraph.graph_colors,
            ref_meta = BuildGraph.ref_meta,

            ref_ids = ref_ids,
            ref_links = BuildLinks.ref_links,
            mm2_indices = BuildLinks.minimap2_index,

            gcs_output_dir = gcs_output_dir
    }
}



