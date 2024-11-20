version 1.0

import "../../structs/Structs.wdl"

task GetRunInfo {

    meta{
        description: "Get ONT run info from a final summary file."
    }

    parameter_meta {
        final_summary: "Sequencing summary file."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        String final_summary

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        gsutil cat "~{final_summary}" | sed 's/=[[:space:]]*$/=unknown/' | sed 's/=/\t/g' > run_info.txt
    >>>

    output {
        Map[String, String] run_info = read_map("run_info.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task ListFiles {

    meta {
        description: "List files in a GCS directory."
    }

    parameter_meta {
        sequencing_summary: "Sequencing summary file."
        suffix: "Suffix of files to list."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        String sequencing_summary
        String suffix

        RuntimeAttr? runtime_attr_override
    }

    String indir = sub(sub(sequencing_summary, basename(sequencing_summary), ""), "/$", "")

    command <<<
        set -euxo pipefail

        gsutil ls "~{indir}/**.~{suffix}*" | grep -v fail > files.txt
        cat files.txt | wc -l > lc.txt
    >>>

    output {
        File manifest = "files.txt"
        Int count = read_int("lc.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task PartitionManifest {

    meta {
        description: "Partition a manifest into chunks."
    }

    parameter_meta {
        manifest: "Manifest to partition."
        N: "Number of chunks to partition into."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        File manifest
        Int N

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        split -a 5 -d --additional-suffix=.txt -e -n l/~{N} ~{manifest} manifest_chunk_
    >>>

    output {
        Array[File] manifest_chunks = glob("manifest_chunk_*.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

