version 1.0

import "Structs.wdl"

# TODO: describe purpose
task MergeBams {
    input {
        Array[File] aligned_shards
        String merged_name
        
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3*ceil(size(aligned_shards, "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx4g -jar /usr/local/bin/gatk.jar MergeSamFiles -I ~{sep=" -I " aligned_shards} -O ~{merged_name} -AS --CREATE_INDEX --USE_JDK_DEFLATER --USE_JDK_INFLATER --VALIDATION_STRINGENCY SILENT
    >>>

    output {
        File merged = "~{merged_name}"
        File merged_bai = basename(merged_name, ".bam") + ".bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2, 
        mem_gb:             20, 
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "kgarimella/lr-align:0.01.17"
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