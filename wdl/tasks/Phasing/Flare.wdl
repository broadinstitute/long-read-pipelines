version 1.0

import "../../structs/Structs.wdl"


task Flare {

    meta {
        description: "Local Ancestry Inference"
    }


    input {
        File ref_vcf
        File ref_vcf_index
        File ref_panel
        File test_vcf
        File test_vcf_index
        File plink_map
        String output_prefix
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        RuntimeAttr? runtime_attr_override
    }


    command <<<
        set -euxo pipefail

        java -jar /LAI/flare.jar ref=~{ref_vcf} gt=~{test_vcf} map=~{plink_map} ref-panel=~{ref_panel} out=~{output_prefix}

    >>>

    output {
        Array[File] output_files = glob("*")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             64,
        disk_gb:            200,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "hangsuunc/flare:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
