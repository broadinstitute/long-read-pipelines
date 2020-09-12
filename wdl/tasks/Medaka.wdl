version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.33/wdl/tasks/Structs.wdl"

task MedakaPolish {
    input {
        File basecalled_reads
        File draft_assembly
        String model #Run `medaka tools list_models` and pick string with the correct pore type, machine, and guppy version

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2 * ceil(size(basecalled_reads, "GB") + size(draft_assembly, "GB"))

    command <<<
        source /medaka/venv/bin/activate

        set -euxo pipefail

        medaka_consensus -i ~{basecalled_reads} -d ~{draft_assembly} -o output -t 8 -m ~{model}
    >>>

    output {
        File polished_assembly = "./output/consensus.fasta"
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:      8,
        mem_gb:         24,
        disk_gb:        disk_size,
        boot_disk_gb:   10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-medaka:0.1.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:        select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:     select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:     select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:        select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:         select_first([runtime_attr.max_retries, default_attr.max_retries])
        gpuType:                "nvidia-tesla-t4"
        gpuCount:               1
        nvidiaDriverVersion:    "418.87.00"
        zones:                  ["us-east1-c"]
        cpuPlatform:            "Intel Haswell"
        docker:             select_first([runtime_attr.docker, default_attr.docker])
    }
}
