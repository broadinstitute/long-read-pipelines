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

    String ref_filename = basename(ref_path)

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
        ln -s ~{ref_path} ~{ref_filename}

        # Build minimap2 index for ref contig mapping in later steps
        minimap2 -x asm5 -d ~{ref_filename}.mm2 ~{ref_path}

        cd ../..

        # Build actual links
        fusilli db build-links fusilli_db/ ~{ref_id}
    >>>

    output {
        File minimap2_index = "fusilli_db/~{ref_id}/~{ref_filename}.mm2"
        File ref_links = "fusilli_db/~{ref_id}/~{ref_filename}.links"
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

        # Symlink individual reference files
        while IFS= read -r line; do
            parts=(${line})

            mkdir -p "${parts[0]}"
            ln -s "${parts[1]}" "${parts[0]}/${parts[1]##*/}"
            ln -s "${parts[2]}" "${parts[0]}/${parts[2]##*/}"
        done < <(paste ~{write_lines(ref_ids)} ~{write_lines(ref_links)} ~{write_lines(mm2_indices)})

        gsutil -m cp -r . ~{output_dir}
    >>>

    output {
        String ref_graph = "~{output_dir}/reference_graph.gfa"
        String ref_graph_colors = "~{output_dir}/reference_graph.bfg_colors"
        String ref_meta = "~{output_dir}/reference_meta.tsv"
    }

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


task BuildSampleGraph {
    input {
        String sample_id

        File reads_fq1
        File? reads_fq2

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        sample_id: "Name of the sample"
        reads_fq1: "Path to FASTQ file with reads"
        reads_fq2: "In case of Illumina paired-end reads, the path to the second FASTQ file"
    }

    Array[File] reads = select_all([reads_fq1, reads_fq2])
    Int num_threads = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    command <<<
        # Activate fusilli conda env
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        fusilli sample build-and-clean -t ~{num_threads} ~{sample_id} ~{sep=' ' reads}
    >>>

    output {
        File sample_graph = "~{sample_id}/~{sample_id}.gfa"
        File sample_graph_colors = "~{sample_id}/~{sample_id}.bfg_colors"
        File kmer_counts = "~{sample_id}/~{sample_id}.counts"
        File cleaned_graph = "~{sample_id}/~{sample_id}.cleaned.gfa"
        File cleaned_graph_colors = "~{sample_id}/~{sample_id}.cleaned.bfg_colors"
        File clean_meta = "~{sample_id}/sample_clean_meta.json"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            20,
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
        cpuPlatform:            "Intel Cascade Lake"
    }
}

task ConstructSampleLinks {
    input {
        String sample_id
        File sample_graph
        File sample_graph_colors

        File ref_meta
        Array[String] ref_ids
        Array[File] ref_paths

        File reads_fq1
        File? reads_fq2

        Int prune_threshold

        RuntimeAttr? runtime_attr_override
    }

    Array[File] reads = select_all([reads_fq1, reads_fq2])

    command <<<
        # Activate fusilli conda env
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        # First rebuild DB folder structure
        mkdir -p fusilli_db
        cd fusilli_db

        ln -s ~{ref_meta} reference_meta.tsv

        # Symlink individual reference files
        while IFS= read -r line; do
            parts=(${line})

            mkdir -p "${parts[0]}"
            ln -s "${parts[1]}" "${parts[0]}/${parts[1]##*/}"
        done < <(paste ~{write_lines(ref_ids)} ~{write_lines(ref_paths)})

        cd ..

        # Organize sample related files
        mkdir ~{sample_id}
        cd ~{sample_id}

        ln -s ~{sample_graph} ~{sample_id}.cleaned.gfa
        ln -s ~{sample_graph_colors} ~{sample_id}.cleaned.bfg_colors

        cd ..

        fusilli sample construct-links ~{sample_id} -o ~{sample_id}.links -r fusilli_db -s ~{sep=' ' reads}
        fusilli sample prune-links ~{sample_id}.links -t ~{prune_threshold} --in-place
    >>>

    output {
        File sample_links = "~{sample_id}.links"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            20,
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
        cpuPlatform:            "Intel Cascade Lake"
    }
}

task BuildCombinedGraph {
    input {
        Array[String] sample_ids
        Array[File] cleaned_graphs
        Array[File] cleaned_graph_colors

        RuntimeAttr? runtime_attr_override
    }

    Int num_threads = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    command <<<
        # Activate fusilli conda env
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        # Organize sample graphs
        while IFS= read -r line; do
            parts=(${line})

            mkdir -p "${parts[0]}"
            ln -s "${parts[1]}" "${parts[0]}/${parts[1]##*/}"
            ln -s "${parts[1]}" "${parts[0]}/${parts[2]##*/}"
        done < <(paste ~{write_lines(sample_ids)} ~{write_lines(cleaned_graphs)} ~{write_lines(cleaned_graph_colors)})

        fusilli samples combine **/*.cleaned.gfa -o combined_graph.gfa -t ~{num_threads}
    >>>

    output {
        File combined_graph = "combined_graph.gfa"
        File combined_graph_colors = "combined_graph.bfg_colors"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            50,
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
        cpuPlatform:            "Intel Cascade Lake"
    }
}

task FindVariantKmers {
    input {
        File ref_graph
        File ref_graph_colors

        File combined_graph
        File combined_graph_colors

        RuntimeAttr? runtime_attr_override
    }

    Int num_threads = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    command <<<
        # Activate fusilli conda env
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        fusilli samples scan ~{ref_graph} ~{combined_graph} -o kmers
    >>>

    output {
        File variant_kmers = "kmers.variants.txt"
        File nonref_kmers = "kmers.nonref.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            50,
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
        cpuPlatform:            "Intel Cascade Lake"
    }
}

task AssembleVariantContigs {
    input {
        String sample_id
        File cleaned_graph
        File cleaned_graph_colors

        File linkdb

        File variant_kmers
        File nonref_kmers

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        # Activate fusilli conda env
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        mkdir ~{sample_id}
        cd ~{sample_id}

        ln -s ~{cleaned_graph} ~{sample_id}.cleaned.gfa
        ln -s ~{cleaned_graph_colors} ~{sample_id}.cleaned.bfg_colors

        cd ..

        fusilli sample assemble ~{sample_id} ~{variant_kmers} ~{nonref_kmers} -l ~{linkdb} -o ~{sample_id}.contigs.fasta
    >>>

    output {
        File variant_contigs = "~{sample_id}.contigs.fasta"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            20,
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
        cpuPlatform:            "Intel Cascade Lake"
    }
}

task FinalizeAssembly {
    input {
        Array[String] sample_ids
        Array[File] linkdbs
        Array[File] variant_contigs

        File combined_graph
        File combined_graph_colors
        File variant_kmers
        File nonref_kmers

        String gcs_output_dir

        RuntimeAttr? runtime_attr_override
    }

    String output_dir = sub(gcs_output_dir, "/$", "")

    command <<<
        set -euxo pipefail

        # Create final folder structure
        mkdir output
        cd output

        # Organize sample data
        while IFS= read -r line; do
            parts=(${line})

            mkdir -p "${parts[0]}"
            ln -s "${parts[1]}" "${parts[0]}/${parts[1]##*/}"
            ln -s "${parts[1]}" "${parts[0]}/${parts[2]##*/}"
        done < <(paste ~{write_lines(sample_ids)} ~{write_lines(linkdbs)} ~{write_lines(variant_contigs)})

        ln -s ~{combined_graph} combined_graph.gfa
        ln -s ~{combined_graph_colors} combined_graph.bfg_colors
        ln -s ~{variant_kmers} kmers.variants.txt
        ln -s ~{nonref_kmers} kmers.nonref.txt

        gsutil -m cp -r . ~{output_dir}
    >>>

    output {
        String combined_graph = "~{output_dir}/combined_graph.gfa"
        String combined_graph_colors = "~{output_dir}/combined_graph.bfg_colors"
        String variant_kmers = "~{output_dir}/kmers.variants.txt"
        String nonref_kmers = "~{output_dir}/kmers.nonref.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            20,
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
