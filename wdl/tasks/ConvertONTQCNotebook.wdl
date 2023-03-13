version 1.0

import "Structs.wdl"

task ConvertONTQCNotebook {
    input {
        String workspace_namespace
        String workspace_name

        String sample_name
        String notebook_path
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us-central1-docker.pkg.dev/broad-dsp-lrma/general/lr-ontqc:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euxo pipefail

        export WORKSPACE_NAMESPACE="~{workspace_namespace}"
        export WORKSPACE_NAME="~{workspace_name}"

        papermill --parameters SAMPLE "~{sample_name}" ~{notebook_path} "~{sample_name}.qc.ipynb"
        jupyter nbconvert --to html --no-input --no-prompt "~{sample_name}.qc.ipynb"
    >>>

    output {
        File ont_qc_report = "~{sample_name}.qc.html"
    }

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

}
