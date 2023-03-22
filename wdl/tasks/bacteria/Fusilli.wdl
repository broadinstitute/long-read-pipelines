version 1.0

import "../Structs.wdl"


task CreateDatabaseFromPanaroo {
    input {
        String db_name
        File panaroo_graph
        File panaroo_gene_data

        Int index_k = 23
        Int? index_g

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us-central1-docker.pkg.dev/broad-dsp-lrma/fusilli/fusilli:devel"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int num_cpus = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    runtime {
        cpu:                    num_cpus
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

    command <<<
        set -euxo pipefail

        mkdir panaroo
        cd panaroo
        ln -s ~{panaroo_graph} final_graph.gml
        ln -s ~{panaroo_gene_data} gene_data.csv
        cd ..

        fusilli db create -o "~{db_name}" --from-panaroo panaroo
        fusilli db index -j~{num_cpus} -k~{index_k} ~{"-g" + index_g} "~{db_name}"

        tar -caf "~{db_name}.tar.gz" "~{db_name}"
    >>>

    output {
        File fusilli_db_tar = "~{db_name}.tar.gz"
        File fusilli_db_index = "~{db_name}.fidx"
    }
}

task PrepareSample {
    input {
        String sample_name
        File illumina_fq1
        File illumina_fq2

        Int k = 49
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us-central1-docker.pkg.dev/broad-dsp-lrma/fusilli/fusilli:devel"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int num_cpus = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    runtime {
        cpu:                    num_cpus
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

    command <<<
        set -euxo pipefail

        fusilli sample prepare -o "~{sample_name}" -j~{num_cpus} -k~{k} ~{illumina_fq1} ~{illumina_fq2}

        tar -caf "~{sample_name}.tar.gz" -C "~{sample_name}" .
    >>>

    output {
        File fusilli_prepared_sample = "~{sample_name}.tar.gz"
    }
}

task CallGenes {
    input {
        File fusilli_db_tar
        File fusilli_db_index

        File fusilli_prepared_sample

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us-central1-docker.pkg.dev/broad-dsp-lrma/fusilli/fusilli:devel"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int num_cpus = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    runtime {
        cpu:                    num_cpus
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

    command <<<
        set -euxo pipefail

        sample_name="~{basename(fusilli_prepared_sample, ".tar.gz")}"
        mkdir "${sample_name}"
        tar -xaf "~{fusilli_prepared_sample}" -C "${sample_name}"

        tar -xaf ~{fusilli_db_tar}

        fusilli call genes -j~{num_cpus} "~{basename(fusilli_db_tar, ".tar.gz")}" "${sample_name}"

        mv "${sample_name}/gene_families/gene_predictions.tsv" .
    >>>

    output {
        File gene_predictions = "gene_predictions.tsv"
    }
}
