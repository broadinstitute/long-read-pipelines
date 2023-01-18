version 1.0

import "../Structs.wdl"

task BaktaDBDownload {
    input {
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        bakta_db download --output .

        tar -caf baktadb.tar.bz2 -C db .
    >>>

    output {
        File bakta_db = "baktadb.tar.bz2"
    }


    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            200,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "quay.io/biocontainers/bakta:1.6.1--pyhdfd78af_0"
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


task BaktaAnnotate {
    input {
        File bakta_db_tar
        File genome_fasta
        String? fname_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int num_cores = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    String prefix = select_first([fname_prefix, basename(genome_fasta)])

    command <<<
        set -euxo pipefail

        mkdir bakta_db
        tar -xaf ~{bakta_db_tar} -C bakta_db

        mkdir output
        bakta --db bakta_db --output output --complete --threads ~{num_cores} \
            --keep-contig-headers --prefix ~{prefix} --verbose \
            ~{genome_fasta}
    >>>

    output {
        File tsv = "output/~{prefix}.tsv"
        File gff = "output/~{prefix}.gff3"
        File genbank = "output/~{prefix}.gbff"
        File embl = "output/~{prefix}.embl"
        File ffn = "output/~{prefix}.ffn"
        File faa = "output/~{prefix}.faa"
        File hypotheticals_tsv = "output/~{prefix}.hypotheticals.tsv"
        File hypotheticals_faa = "output/~{prefix}.hypotheticals.faa"

        File summary = "output/~{prefix}.txt"
        File plot_png = "output/~{prefix}.png"
        File plot_svg = "output/~{prefix}.svg"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "quay.io/biocontainers/bakta:1.6.1--pyhdfd78af_0"
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
