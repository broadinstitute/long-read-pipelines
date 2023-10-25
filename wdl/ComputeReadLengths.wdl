version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/PBUtils.wdl" as PB
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/Finalize.wdl" as FF

workflow ComputeReadLengths {
    input {
        File input_pbi
        String movie_name

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ComputeReadLengths/"

    call GetReadLengths { input: input_pbi = input_pbi, movie_name = movie_name }

    ##########
    # Finalize
    ##########

    call FF.FinalizeToFile { input: file = GetReadLengths.read_lengths, outdir = outdir }

    output {
        File read_lengths = GetReadLengths.read_lengths
    }
}

task GetReadLengths {
    input {
        File input_pbi
        String movie_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(input_pbi, "GB"))

    command <<<
        set -euxo pipefail

        python3 /usr/local/bin/compute_read_length_hist.py ~{input_pbi} | gzip > ~{movie_name}.read_lengths.txt.gz
    >>>

    output {
        File read_lengths = "~{movie_name}.read_lengths.txt.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.40-kvg"
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
