version 1.0

import "tasks/Structs.wdl" as Structs
import "tasks/Utils.wdl" as Utils

workflow CallTelomeres {
    input {
        String gcs_fast5_dir
    }

    call Utils.ListFilesOfType { input: gcs_dir = gcs_fast5_dir, suffixes = [ ".fast5" ] }
    call Utils.ChunkManifest { input: manifest = ListFilesOfType.manifest, manifest_lines_per_chunk = 1000 }

    scatter (manifest_chunk in ChunkManifest.manifest_chunks) {
        call Bonito {
            input:
                fast5_files = read_lines(manifest_chunk),
        }
    }

    call Utils.MergeFastqGzs { input: fastq_gzs = Bonito.basecalls, prefix = "basecalls" }

    output {
        File basecalls = MergeFastqGzs.merged_fastq_gz
    }
}

task Bonito {
    input {
        Array[File] fast5_files

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(fast5_files, "GB"))

    command <<<
        set -euxo pipefail

        num_cores=$(grep -c '^processor' /proc/cpuinfo | awk '{ print $1 - 1 }')
        fast5_dir=$(dirname ~{fast5_files[0]})

        pwd
        ls

        source /bonito-0.3.0/venv3/bin/activate

        bonito basecaller $model_dir $fast5_dir | gzip > recalled.fa.gz
    >>>

    output {
        File basecalls = "recalled.fa.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             50,
        disk_gb:            disk_size,
        boot_disk_gb:       30,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/ymostovoy/lr-bonito:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        gpuType:                "nvidia-tesla-p100"
        gpuCount:               1
        nvidiaDriverVersion:    "418.152.00"
        zones:                  ["us-central1-c", "us-central1-f", "us-east1-b", "us-east1-c", "us-west1-a", "us-west1-b"]
        cpuPlatform:            "Intel Haswell"
    }
}
