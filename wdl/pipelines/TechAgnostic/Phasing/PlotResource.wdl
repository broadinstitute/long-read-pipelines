version 1.0

import "../../../structs/Structs.wdl"

workflow Statistics {
    meta{
        description : "..."
    }
    parameter_meta {
        
    }

    input {
        File resource_log
    }
    
    call resource_usage { input:
        resource_log = resource_log
    }

    output{
        File resource_plot = resource_usage.resource_plot
    }
}

task resource_usage {

    input {
        File resource_log
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        source("./opt/plot.resources.R")
        resources = get_time_series_df(~{resource_log})
        plot_time_series_df(resources, "resource_plot.pdf", "resources usage by task")

    >>>

    output {
        File resource_plot = "resource_plot.pdf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             64,
        disk_gb:            100,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-resource-visual:latest"
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
