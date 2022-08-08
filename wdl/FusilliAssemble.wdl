version 1.0

import "tasks/Fusilli.wdl" as Fusilli

workflow FusilliAssemble {
    input {
        String fusilli_db_gcs
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

        String gcs_output_dir
    }

    String output_dir = sub(gcs_output_dir, "/$", "")

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

        call Fusilli.AssembleVariantContigs as AssembleVariantContigs {
            input:
                sample_id = sample_ids[i],
                cleaned_graph = cleaned_sample_graphs[i],
                cleaned_graph_colors = cleaned_sample_graph_colors[i],
                linkdb = ConstructSampleLinks.sample_links,
                variant_kmers = FindVariantKmers.variant_kmers,
                nonref_kmers = FindVariantKmers.nonref_kmers
        }

        call Fusilli.BuildRefPanels as BuildRefPanels {
            input:
                ref_db_gcs = fusilli_db_gcs,

                variant_contigs = AssembleVariantContigs.variant_contigs
        }

    }

    call Fusilli.FinalizeAssembly as FinalizeAssembly {
        input:
            gcs_output_dir = output_dir,
            sample_ids = sample_ids,
            linkdbs = ConstructSampleLinks.sample_links,
            variant_contigs = AssembleVariantContigs.variant_contigs,
            combined_graph = BuildCombinedGraph.combined_graph,
            combined_graph_colors = BuildCombinedGraph.combined_graph_colors,
            variant_kmers = FindVariantKmers.variant_kmers,
            nonref_kmers = FindVariantKmers.nonref_kmers
    }

    # Workaround to define finalized paths for sample data
    # WDL doesn't provide functionality to apply a function on an Array (e.g., applying `basename` on all
    # entries in an array), and scatter can't be used in a task, so we do it here.
    scatter(i in range(length(sample_ids))) {
        String finalized_linkdbs = "~{output_dir}/~{sample_ids[i]}/~{sample_ids[i]}.links"
        String finalized_contigs = "~{output_dir}/~{sample_ids[i]}/~{sample_ids[i]}.contigs.fasta"

        call Fusilli.FinalizeRefPanels {
            input:
                sample_id = sample_ids[i],
                ref_panels = BuildRefPanels.ref_panels[i],

                gcs_output_dir = output_dir
        }

        String finalized_ref_panels = "~{output_dir}/~{sample_ids[i]}/ref_panels"
    }

    output {
        String combined_graph = FinalizeAssembly.combined_graph
        String combined_graph_colors = FinalizeAssembly.combined_graph_colors
        String variant_kmers = FinalizeAssembly.variant_kmers
        String nonref_kmers = FinalizeAssembly.nonref_kmers
        Array[String] linkdbs = finalized_linkdbs
        Array[String] variant_contigs = finalized_contigs
        Array[String] ref_panels = finalized_ref_panels
    }
}
