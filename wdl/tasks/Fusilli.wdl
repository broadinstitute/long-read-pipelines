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

    Int disk_size = 1 + ceil((length(ref_ids) * 25) / 1024)

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
        fusilli sample kmer-spectrum ~{sample_id}/~{sample_id}.counts -o ~{sample_id}/~{sample_id}.spectrum.png
    >>>

    output {
        File sample_graph = "~{sample_id}/~{sample_id}.gfa"
        File sample_graph_colors = "~{sample_id}/~{sample_id}.bfg_colors"
        File kmer_counts = "~{sample_id}/~{sample_id}.counts"
        File kmer_spectrum = "~{sample_id}/~{sample_id}.spectrum.png"
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
    }
}

task ConstructSampleLinks {
    input {
        String sample_id
        File sample_graph
        File sample_graph_colors

        File reads_fq1
        File? reads_fq2

        RuntimeAttr? runtime_attr_override
    }

    Array[File] reads = select_all([reads_fq1, reads_fq2])

    command <<<
        # Activate fusilli conda env
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        # Organize sample related files
        mkdir ~{sample_id}
        cd ~{sample_id}

        ln -s ~{sample_graph} ~{sample_id}.cleaned.gfa
        ln -s ~{sample_graph_colors} ~{sample_id}.cleaned.bfg_colors

        cd ..

        fusilli sample links-from-reads ~{sample_id} ~{sep=' ' reads} -o ~{sample_id}.links \
            --prune-first-pass-db 3
    >>>

    output {
        File sample_links = "~{sample_id}.links"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             32,
        disk_gb:            20,
        boot_disk_gb:       10,
        preemptible_tries:  1,
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
    }
}

task MergeSampleAndRefLinks {
    input {
        String sample_id
        File cleaned_graph
        File cleaned_graph_colors
        File sample_links

        File ref_meta
        Array[String] ref_ids
        Array[File] ref_fastas

        Int prune_threshold

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        # Activate fusilli conda env
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        mkdir -p fusilli_db
        cd fusilli_db

        ln -s ~{ref_meta} reference_meta.tsv

        # Symlink individual reference files
        while IFS=$'\t' read -r ref_id ref_fasta; do
            mkdir -p "${ref_id}"
            ln -s "${ref_fasta}" "${ref_id}/${ref_fasta##*/}"
        done < <(paste ~{write_lines(ref_ids)} ~{write_lines(ref_fastas)})

        cd ..

        mkdir ~{sample_id}
        cd ~{sample_id}

        ln -s ~{cleaned_graph} ~{sample_id}.cleaned.gfa
        ln -s ~{cleaned_graph_colors} ~{sample_id}.cleaned.bfg_colors
        ln -s ~{sample_links} ~{sample_id}.links

        cd ..

        fusilli sample prune-links ~{sample_id}/~{sample_id}.links -t ~{prune_threshold} -o ~{sample_id}.pruned.links

        fusilli sample links-from-refs ~{sample_id} fusilli_db/ -o ~{sample_id}.refs.links
        fusilli sample merge-links ~{sample_id}.refs.links ~{sample_id}.pruned.links -o ~{sample_id}.temp.links
    >>>

    output {
        File sample_ref_links = "~{sample_id}.refs.links"
        File sample_merged_links = "~{sample_id}.links"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
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

        fusilli sample assemble ~{sample_id} ~{variant_kmers} ~{nonref_kmers} -l ~{linkdb} -o ~{sample_id}.tmp.fasta
        fusilli sample rm-contained ~{sample_id}.tmp.fasta -o ~{sample_id}.contigs.fasta > num_contained
    >>>

    output {
        File variant_contigs = "~{sample_id}.contigs.fasta"
        Int num_contained = read_int("num_contained")
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
    }
}

task BuildRefPanels {
    input {
        String ref_db_gcs

        File variant_contigs

        RuntimeAttr? runtime_attr_override
    }

    String db_name = basename(ref_db_gcs)

    command <<<
        # Activate fusilli conda env
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        mkdir fusilli_db
        gsutil -m cp -r ~{ref_db_gcs}/ .

        fusilli call build-ref-panels ~{db_name} ~{variant_contigs} -o ref_contigs/
    >>>

    output {
        Array[File] ref_panels = glob("ref_contigs/*.ref_contigs.fa")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
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

    String output_dir = sub(gcs_output_dir, "/$", "") + "/assembly"

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
            ln -s "${parts[2]}" "${parts[0]}/${parts[2]##*/}"
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

task FinalizeRefPanels {
    input {
        String sample_id
        Array[File] ref_panels

        String gcs_output_dir

        RuntimeAttr? runtime_attr_override
    }

    String output_dir = sub(gcs_output_dir, "/$", "") + "/assembly"
    String sample_gcs_dir = output_dir + "/" + sample_id

    command <<<
        set -euxo pipefail

        mkdir ref_panels/
        cd ref_panels

        while IFS= read -r line; do
            fname=${line##*/}
            ln -s ${line} ${fname}
        done < ~{write_lines(ref_panels)}

        cd ..

        gsutil -m cp -r ref_panels ~{sample_gcs_dir}
    >>>

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

task ChunkSampleContigs {
    input {
        File sample_contigs

        Int? chunk_size = 1000
        RuntimeAttr? runtime_attr_override
    }

    String contig_basename = basename(sample_contigs)

    command <<<
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        ln -s ~{sample_contigs} ~{contig_basename}

        fusilli utils split-fasta ~{contig_basename} -c ~{chunk_size}
    >>>

    output {
        Array[File] chunks = glob("*.chunk*")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
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

task TesseraeAlign {
    input {
        String sample_id
        File sample_contigs
        String fusilli_run_gcs
        File? hmm_config

        RuntimeAttr? runtime_attr_override
    }

    String sample_gcs_dir = fusilli_run_gcs + "/" + sample_id

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        3,
        docker:             "us-east1-docker.pkg.dev/broad-dsp-lrma/fusilli/tesserae2:devel"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int max_memory_tesserae = round(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 2)

    command <<<
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        gsutil -m cp -r ~{sample_gcs_dir} .

        mkdir output
        tesserae multi-align "~{sample_id}/ref_panels/{query_id}.ref_contigs.fa" ~{sample_contigs} \
            ~{if defined(hmm_config) then "-c ~{hmm_config}" else ""} -O bam -o output/ -m ~{max_memory_tesserae}G
    >>>

    output {
        Array[File] alignments = glob("output/*.bam")
        File unaligned = "output/unaligned.fasta.gz"
    }

    #########################
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

task MakeCalls {
    input {
        String fusilli_db_gcs
        Array[String] ref_ids
        Array[File] ref_fastas

        String fusilli_run_gcs
        String sample_id
        Array[File] aligned_sample_contigs

        RuntimeAttr? runtime_attr_override
    }

    String db_name = basename(fusilli_db_gcs)
    String sample_run_data = fusilli_run_gcs + "/assembly/" + sample_id

    command <<<
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        mkdir db
        mkdir run
        gsutil -m cp -r ~{fusilli_db_gcs} db/
        gsutil -m cp -r ~{sample_run_data} run/

        # Create symlinks to actual reference FASTAs
        cd db/~{db_name}

        while IFS=$'\t' read -r ref_id ref_fasta; do
            ln -s "${ref_fasta}" "${ref_id}/${ref_fasta##*/}"
        done < <(paste ~{write_lines(ref_ids)} ~{write_lines(ref_fastas)})

        cd ../..

        mkdir output
        fusilli call run db/~{db_name} ~{sample_id} run/~{sample_id}/ref_panels \
            -o output/ < ~{write_lines(aligned_sample_contigs)}

        # Sort, normalize and deduplicate VCF. We output to BCF to be able to get indexed by bcftools index in
        # later Fusilli stage
        for vcf in output/*.vcf; do
            filename=${vcf##*/}
            ref_name=${filename%.vcf}
            bcftools sort -Ou ${vcf} \
                | bcftools norm --fasta-ref $(fusilli db get-ref-path db/~{db_name} ${ref_name}) -d exact \
                    --check-ref wx -Ob -o ${vcf%.vcf}.bcf
        done
    >>>

    output {
        File alignment_stats = "output/~{sample_id}.aln_stats.tsv"
        Array[File] aligned_ref_contigs = glob("output/*.bam")
        Array[File] aligned_ref_contigs_bai = glob("output/*.bai")
        Array[File] per_ref_vcfs = glob("output/*.bcf")
        Array[File] per_ref_bed = glob("output/*.bed")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            30,
        boot_disk_gb:       10,
        preemptible_tries:  2,
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

task FinalizeCalls {
    input {
        String fusilli_run_gcs
        String call_id
        String sample_id

        File alignment_stats
        Array[File] aligned_sample_contigs
        Array[File] unaligned_fastas
        Array[File] aligned_ref_contigs
        Array[File] aligned_ref_contigs_bai
        Array[File] per_ref_vcfs
        Array[File] per_ref_bed

        RuntimeAttr? runtime_attr_override
    }

    String call_output_dir = fusilli_run_gcs + "/calls/" + call_id

    command <<<
        mkdir ~{sample_id}
        cd ~{sample_id}

        ln -s ~{alignment_stats} ~{basename(alignment_stats)}

        mkdir aligned_sample_contigs
        mkdir aligned_ref_contigs
        mkdir per_ref_vcfs
        mkdir per_ref_bed

        cat $(< ~{write_lines(unaligned_fastas)}) > unaligned.fasta.gz

        for sample_contig in $(< ~{write_lines(aligned_sample_contigs)}); do
            ln -s "${sample_contig}" "aligned_sample_contigs/${sample_contig##*/}"
        done

        while IFS=$'\t' read -r ref_contig_bam ref_contig_bai; do
            ln -s "${ref_contig_bam}" "aligned_ref_contigs/${ref_contig_bam##*/}"
            ln -s "${ref_contig_bai}" "aligned_ref_contigs/${ref_contig_bai##*/}"
        done < <(paste ~{write_lines(aligned_ref_contigs)} ~{write_lines(aligned_ref_contigs_bai)})

        for vcf in $(< ~{write_lines(per_ref_vcfs)}); do
            ln -s "${vcf}" "per_ref_vcfs/${vcf##*/}"
        done

        for bed in $(< ~{write_lines(per_ref_bed)}); do
            ln -s "${bed}" "per_ref_bed/${bed##*/}"
        done

        # First a non-recursive copy to ensure all folders exists
        gsutil cp ~{basename(alignment_stats)} ~{call_output_dir}/~{sample_id}/~{basename(alignment_stats)}

        cd ..

        gsutil -m cp -r ~{sample_id} ~{call_output_dir}
    >>>

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


task CreatePseudoRef {
    input {
        File ref_meta

        Array[String] ref_ids
        Array[File] ref_fastas
        Array[File] ref_gffs

        Array[File] alignment_stats

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        mkdir db
        cd db

        ln -s ~{ref_meta} reference_meta.tsv

        while IFS=$'\t' read -r ref_id ref_fasta ref_gff; do
            mkdir "${ref_id}"

            ln -s "${ref_fasta}" "${ref_id}/${ref_fasta##*/}"
            ln -s "${ref_gff}" "${ref_id}/${ref_gff##*/}"
        done < <(paste ~{write_lines(ref_ids)} ~{write_lines(ref_fastas)} ~{write_lines(ref_gffs)})

        cd ..

        fusilli call build-pseudo-ref db/ $(< ~{write_lines(alignment_stats)}) -o pseudo_ref.fasta
        samtools faidx pseudo_ref.fasta
    >>>

    output {
        File pseudo_ref = "pseudo_ref.fasta"
        File pseudo_ref_fai = "pseudo_ref.fasta.fai"
        File pseudo_gff = "pseudo_ref.gff"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
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
    }
}

task TranslateCoordinatesToPseudoRef {
    input {
        File pseudo_ref

        String sample_id

        Array[File] per_ref_vcfs
        Array[File] per_ref_bed

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        mkdir to_pseudo_vcfs
        mkdir to_pseudo_beds

        # Translate coordinates to pseudo-ref coordinates
        for vcf in $(< ~{write_lines(per_ref_vcfs)}); do
            vcf_filename=${vcf##*/}
            ref_name=${vcf_filename%.bcf}
            fusilli call translate-to-pseudo-ref ~{pseudo_ref} -m vcf -r ${ref_name} \
                -o to_pseudo_vcfs/${vcf_filename} ${vcf}

            # Index for bcftools concat down below
            bcftools index to_pseudo_vcfs/${vcf_filename}
        done

        for bed in $(< ~{write_lines(per_ref_bed)}); do
            bed_filename=${bed##*/}
            ref_name=${bed_filename%.bed}
            fusilli call translate-to-pseudo-ref ~{pseudo_ref} -m bed -r ${ref_name} \
                -o to_pseudo_beds/${bed_filename} ${bed}
        done

        # Combine all per-ref VCFs to a single VCF for this sample
        bcftools concat to_pseudo_vcfs/*.bcf -Ob -o ~{sample_id}.bcf
        bcftools index ~{sample_id}.bcf

        # Combine all per-ref BEDs to a single BED file
        cat to_pseudo_beds/*.bed > ~{sample_id}.bed
    >>>

    output {
        File pseudo_vcf = "~{sample_id}.bcf"
        File pseudo_vcf_csi = "~{sample_id}.bcf.csi"
        File pseudo_bed = "~{sample_id}.bed"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             16,
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
    }
}

task FinalizePseudoRef {
    input {
        String fusilli_run_gcs
        String call_id

        File pseudo_ref
        File pseudo_ref_fai
        File pseudo_gff

        Array[File] pseudo_vcfs
        Array[File] pseudo_vcf_csis
        Array[File] pseudo_beds

        RuntimeAttr? runtime_attr_override
    }

    String output_dir = fusilli_run_gcs + "/calls/" + call_id + "/pseudo_ref"

    command <<<
        gsutil cp ~{pseudo_ref} ~{pseudo_ref_fai} ~{pseudo_gff} ~{output_dir}
        gsutil -m cp $(< ~{write_lines(pseudo_vcfs)}) ~{output_dir}/sample_data
        gsutil -m cp $(< ~{write_lines(pseudo_vcf_csis)}) ~{output_dir}/sample_data
        gsutil -m cp $(< ~{write_lines(pseudo_beds)}) ~{output_dir}/sample_data
    >>>

    output {
        String pseudo_ref = "~{output_dir}/~{basename(pseudo_ref)}"
        String pseudo_gff = "~{output_dir}/~{basename(pseudo_gff)}"
        String sample_data_dir = "~{output_dir}/sample_data"
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
