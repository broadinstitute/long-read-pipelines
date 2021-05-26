version 1.0

######################################################################################
## A pipeline for running the Winnowmap aligner
######################################################################################

import "Structs.wdl"

task Align {
    input {
        File reads
        String genome_size
        String prefix
        String preset

        Int? min_depth_edge = 3

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:           "Raw reads in either fa or fq format"
        genome_size:     "Estimated genome size, can use k/m/g suffixes (e.g. 3g for the human genome)"
        prefix:          "Output file prefix"
        preset:          "data preset (\"rs\" for PacBio RSII, \"sq\" for PacBio Sequel, \"ccs\" for PacBio CCS reads and \"ont\" for Oxford Nanopore)"

        min_depth_edge:  "Minimum read support for valid edges. Default is 3, can set to 2 for low sequence depth or 4 for very high sequence depth"
    }

    Int disk_size = 10 * ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        wtdbg2 -t $num_core -x ~{preset} -g ~{genome_size} -e ~{min_depth_edge} -i ~{reads} -fo ~{prefix}
        wtpoa-cns -t $num_core -i ~{prefix}.ctg.lay.gz -fo ~{prefix}.ctg.fa
    >>>

    output {
 #       File lay = "~{prefix}.ctg.lay.gz"
        File fa = "~{prefix}.ctg.fa"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-winnowmap:2.03"
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