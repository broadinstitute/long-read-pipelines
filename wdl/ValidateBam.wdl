version 1.0

import "Structs.wdl"

# TODO: describe purpose
task ValidateBam {
    input {
        File input_bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(1.2*size(input_bam, "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx4g -jar /usr/local/bin/gatk.jar ValidateSamFile -I ~{input_bam} -O bam_validation_report.txt --IGNORE_WARNINGS
    >>>

    output {
        File report = "bam_validation_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2, 
        mem_gb:             8, 
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