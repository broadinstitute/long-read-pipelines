version 1.0

import "Structs.wdl"

task RepeatMasker {
    input {
        File fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 5*ceil(size(fasta, "GB"))+20

    String prefix = basename(fasta)

    command <<<
        set -euxo pipefail
  
        RepeatMasker -e rmblast -pa 4 -s -species human ~{fasta}
        mv ~{fasta}.out ./~{prefix}.out
     >>>

     output {
         File RMout = '~{prefix}.out'
     }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "dfam/tetools:1.8"
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
