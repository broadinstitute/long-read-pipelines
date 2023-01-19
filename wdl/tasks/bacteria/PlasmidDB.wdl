version 1.0

import "../Structs.wdl"


task UpdatePlasmidDB {
    input {
        String plsdb_fasta_url
        String plsdb_meta_url

        String gcs_output_directory

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        curl -L ~{plsdb_fasta_url} -o plsdb.fna.bz2
        curl -L ~{plsdb_meta_url} -o plsdb_meta.tar.bz2

        mkdir plsdb
        tar -xjf plsdb_meta.tar.bz2 -C plsdb/

        mkdir plsdb/plasmids
        split_plsdb_seq.py plsdb.fna.bz2 plsdb/plasmids/
        prepare_plsdb_terra.py plsdb/plsdb.tsv ~{gcs_output_directory}/plasmids/ > plsdb/plsdb.terra.tsv

        gsutil -m rsync -r plsdb ~{gcs_output_directory}
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            20,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us-east1-docker.pkg.dev/broad-dsp-lrma/fusilli/plsdb:latest"
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
