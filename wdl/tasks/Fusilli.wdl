version 1.0

import "Structs.wdl"

task BuildGraph {
    input {
        Array[String] ref_ids
        Array[File] ref_paths

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        ref_ids: "List of reference identifiers to include"
        ref_paths: "List of paths to the reference FASTA files"
    }

    command <<<
        # Activate fusilli conda env
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        echo "ID" > ref_ids.txt
        echo "fpath" > ref_paths.txt

        cat ~{write_lines(ref_ids)} >> ref_ids.txt
        cat ~{write_lines(ref_paths)} >> ref_paths.txt

        paste ref_ids.txt ref_paths.txt \
            | fusilli db create -t ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} --symlink fusilli_db/
    >>>

    output {
        File graph = "fusilli_db/reference_graph.gfa"
        File graph_colors = "fusilli_db/reference_graph.bfg_colors"
        File ref_meta = "fusilli_db/reference_meta.tsv"
    }

    Int disk_size = 1 + ceil((length(ref_ids) * 10) / 1024)

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us-east1-docker.pkg.dev/broad-dsp-lrma/fusilli/fusilli:devel"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task BuildRefLinks {
    input {
        File graph
        File graph_colors
        File ref_meta

        String ref_id
        File ref_path

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        graph: "Bifrost graph GFA"
        graph_colors: "Bifrost graph colors file"
        ref_meta: "Fusilli reference meta sheet"

        ref_id: "Reference ID to build links for"
        ref_path: "Reference FASTA path"
    }

    command <<<
        # Activate fusilli conda env
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        # Mimic Fusilli DB folder structure
        mkdir -p fusilli_db
        cd fusilli_db
        ln -s ~{graph} reference_graph.gfa
        ln -s ~{graph_colors} reference_graph.bfg_colors
        ln -s ~{ref_meta} reference_meta.tsv

        mkdir -p ~{ref_id}
        cd ~{ref_id}
        ln -s ~{ref_path} ~{ref_id}.fa.gz

        # Build minimap2 index for ref contig mapping in later steps
        minimap2 -x asm5 -d ~{ref_id}.fa.gz.mm2 ~{ref_path}

        cd ../..

        # Build actual links
        fusilli db build-links fusilli_db/ ~{ref_id}
    >>>

    output {
        File minimap2_index = "fusilli_db/~{ref_id}/~{ref_id}.fa.gz.mm2"
        File ref_links = "fusilli_db/~{ref_id}/~{ref_id}.links"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             32,
        disk_gb:            5,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us-east1-docker.pkg.dev/broad-dsp-lrma/fusilli/fusilli:devel"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


task FinalizeDB {
    input {
        File graph
        File graph_colors
        File ref_meta

        Array[String] ref_ids
        Array[File] ref_links
        Array[File] mm2_indices

        String gcs_output_dir

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        graph: "Reference graph GFA"
        graph_colors: "Bifrost colors file"
        ref_meta: "Fusilli database meta TSV"

        ref_ids: "List of reference IDs"
        ref_links: "List of paths to the reference Link databases"
        mm2_indices: "List of minimap2 indices for each reference"

        gcs_output_dir: "Final destination directory on GCS to store the Fusilli DB."
    }

    String output_dir = sub(gcs_output_dir, "/$", "")

    command <<<
        set -euxo pipefail

        # Mimic Fusilli DB folder structure
        mkdir -p fusilli_db
        cd fusilli_db

        ln -s ~{graph} reference_graph.gfa
        ln -s ~{graph_colors} reference_graph.bfg_colors
        ln -s ~{ref_meta} reference_meta.tsv

        # Add individual reference files
        while IFS= read -r line; do
            parts=(${line})

            mkdir -p "${parts[0]}"
            ln -s "${parts[1]}" "${parts[0]}/${parts[0]}.links"
            ln -s "${parts[2]}" "${parts[0]}/${parts[0]}.fa.gz.mm2"
        done < <(paste ~{write_lines(ref_ids)} ~{write_lines(ref_links)} ~{write_lines(mm2_indices)})

        cd ..

        gsutil -m cp fusilli_db ~{output_dir}
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
