version 1.0

import "../../../structs/Structs.wdl"
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow BackupWorkspace {
    meta {
        author: "Jonn Smith"
        description: "Backup a Terra workspace, copying the backup to the workspace itself (in a backups folder), and optionally to another folder."
    }

    parameter_meta {
        additional_output_path:   "GCP folder in which to place an additional copy of the backup."
    }

    input {
        String? additional_output_path
    }

    # Get our workspace info:
    call GetWorkspaceInfo as t001_GetWorkspaceInfo { input: }

    # Backup our workspace:
    call RunBackup as t002_RunBackup {
        input:
            namespace = t001_GetWorkspaceInfo.namespace,
            workspace = t001_GetWorkspaceInfo.workspace,
            default_bucket = t001_GetWorkspaceInfo.default_bucket,
    }

    # Copy to another location as well if we have to:
    if (defined(additional_output_path)) {
        call FF.FinalizeToDir as CopyBackupToAlternateLocation {
            input:
                outdir = select_first([additional_output_path]),
                files = [t002_RunBackup.backup_path]
        }
    }
}

task GetWorkspaceInfo {
    meta {
        description: "Get the baseline info for the current workspace.  Assumes that we're running in Terra."
    }

    input {
        RuntimeAttr? runtime_attr_override
    }

    String namespace_file = "namespace.txt"
    String workspace_file = "workspace.txt"
    String default_bucket_file = "default_bucket.txt"

    command <<<
    echo $WORKSPACE_NAMESPACE > ~{namespace_file}
    echo $WORKSPACE_NAME > ~{workspace_file}
    echo $WORKSPACE_BUCKET > ~{default_bucket_file}
    >>>

    output {
        String namespace = read_string(namespace_file)
        String workspace = read_string(workspace_file)
        String default_bucket = read_string(default_bucket_file)
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "ubuntu:20.04"
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

task RunBackup {

    meta {
        description: "Run the backup notebook on the workspace with the given info."
    }
    parameter_meta {
        namespace: "Namespace to which the workspace belongs."
        workspace: "Name of the workspace to backup."
        default_bucket: "Default bucket of the workspace.  This will be used as an initial output folder."
    }

    input {
        String namespace
        String workspace
        String default_bucket

        RuntimeAttr? runtime_attr_override
    }

    String out_notebook_file_name = "out_notebook_filename.txt"
    String out_notebook_html_file_name = "out_notebook_html_filename.txt"
    String backup_path_file = "backup_path.txt"

    command <<<
        set -euxo pipefail

        # Create output notebook name:
        out_notebook_name=$(date +%Y%m%dT%H%M%S)_backup_executed_notebook.ipynb

        # Save output notebook name so we can refer to it in our outputs:
        echo ${out_notebook_name} > ~{out_notebook_file_name}
        echo ${out_notebook_name} | sed 's@\.ipynb@.html@' > ~{out_notebook_html_file_name}

        papermill /00_backup_workspace.ipynb ${out_notebook_name}

        jupyter nbconvert --to html ${out_notebook_name}

        # Get the backup path from the output:
        grep 'Backup location:' $(echo ${out_notebook_name} | sed 's@\.ipynb@.html@') > ~{backup_path_file}
    >>>

    output {
        File output_notebook = read_string(out_notebook_file_name)
        File output_notebook_html = read_string(out_notebook_html_file_name)

        String backup_path = read_string(backup_path_file)
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             10,
        disk_gb:            20,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-backup-workspace:0.0.1"
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
