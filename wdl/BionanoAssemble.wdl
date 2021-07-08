version 1.0

import "Structs.wdl"

workflow CallAssemble {
    input {
        File ref_cmap
        File input_bnx
        File optArguments

        RuntimeAttr? runtime_attr_override
    }

    call Assemble {
        input:
            ref_cmap = ref_cmap,
            input_bnx = input_bnx,
            optArguments = optArguments
    }
}

task Assemble {
    input {
        File ref_cmap
        File input_bnx
        File optArguments

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 50 + ceil(size([ref_cmap, input_bnx], "GB"))

    command <<<
        set -euxo pipefail

        source activate bionano_minimal
        python /home/bionano2/tools/pipeline/1.0/Pipeline/1.0/pipelineCL.py -l /home/bionano2/output -t /home/bionano2/tools/pipeline/1.0/RefAligner/1.0 \
        -b ~{input_bnx} -a ~{optArguments} -r ~{ref_cmap} \
        -y -d -U -i 5 -F 1 -W 1 -c 0 --docker --autoRestart \
        -T 64 -je 64
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          64,
        mem_gb:             256,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/ymostovoy/bionano-pipeline:latest"
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