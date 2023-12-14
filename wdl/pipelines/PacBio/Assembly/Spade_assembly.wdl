version 1.0

import "../../../structs/Structs.wdl"

workflow IllumAssembleWithSPAdes {
    input {
        File illumina_fq1
        File illumina_fq2
    }

    call SPAdesAssemble {
        input:
            illumina_fq1=illumina_fq1,
            illumina_fq2=illumina_fq2
    }

    output {
        Array[File] asmfiles = SPAdesAssemble.assemblies
    }
}

task SPAdesAssemble {
    input {
        File illumina_fq1
        File? illumina_fq2

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            25,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "hangsuunc/spades:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int max_mem = round(select_first([runtime_attr.mem_gb, default_attr.mem_gb])) - 2  # 2 GB buffer
    Int num_threads = round(select_first([runtime_attr.cpu_cores, default_attr.cpu_cores]))

    command <<<
        set -euxo pipefail

        spades.py -t ~{num_threads} -m ~{max_mem} --isolate -1 ~{illumina_fq1} ~{"-2 " + illumina_fq2} -o output
    >>>

    output {
        Array[File] assemblies = glob("*")
    }

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

