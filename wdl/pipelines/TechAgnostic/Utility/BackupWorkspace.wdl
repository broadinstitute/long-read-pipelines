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
    call RunBackupWithPapermill as t002_RunBackupPapermill {
        input:
            namespace = t001_GetWorkspaceInfo.namespace,
            workspace = t001_GetWorkspaceInfo.workspace,
            default_bucket = t001_GetWorkspaceInfo.default_bucket,
    }

    # Backup our workspace:
    call RunBackupWithPython as t003_RunBackupPython {
        input:
            namespace = t001_GetWorkspaceInfo.namespace,
            workspace = t001_GetWorkspaceInfo.workspace,
            default_bucket = t001_GetWorkspaceInfo.default_bucket,
    }

    # Copy to another location as well if we have to:
    if (defined(additional_output_path)) {
        call FF.FinalizeToDir as CopyNotebookBackupToAlternateLocation {
            input:
                outdir = select_first([additional_output_path]),
                files = [t002_RunBackupPapermill.backup_path]
        }

        call FF.FinalizeToDir as CopyPythonBackupToAlternateLocation {
            input:
                outdir = select_first([additional_output_path]),
                files = [t003_RunBackupPython.backup_path]
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

    echo "Namespace: $WORKSPACE_NAMESPACE"
    echo "Workspace: $WORKSPACE_NAME"
    echo "Default Bucket: $WORKSPACE_BUCKET"
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

task RunBackupWithPapermill {

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


task RunBackupWithPython {

    meta {
        description: "Backup the workspace with the given info."
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

    String backup_path_file = "backup_path.txt"

    command <<<
        set -euxo pipefail

        python3 <<CODE

        ################################################################################

        import time
        import os
        import datetime
        import gzip
        import io
        import json

        import pandas as pd
        import firecloud.api as fapi
        import numpy as np

        from google.cloud import storage

        import lrmaCU.terra.table_utils as lrma_table_utils
        import lrmaCU.terra.workspace_utils as lrma_workspace_utils

        ################################################################################

        # Track start time so we can calculate elapsed time for backup:
        START_TIME = time.time()

        # Get the workspace info:
        namespace = os.environ['WORKSPACE_NAMESPACE']
        workspace = os.environ['WORKSPACE_NAME']
        default_bucket = os.environ['WORKSPACE_BUCKET']

        print(f"Namespace: {namespace}")
        print(f"Workspace: {workspace}")
        print(f"Default Bucket: {default_bucket}")

        ################################################################################

        def load_table(namespace, workspace, table_name, store_membership=False):
            ent_old = fapi.get_entities(namespace, workspace, table_name).json()
            tbl_old = None

            membership = None
            if len(ent_old) > 0:
                tbl_old = pd.DataFrame(list(map(lambda e: e['attributes'], ent_old)))
                tbl_old[f"entity:{table_name}_id"] = list(map(lambda f: f['name'], ent_old))

                if store_membership:
                    membership = list(map(lambda g: set(map(lambda h: h['entityName'], g['items'])), tbl_old['samples']))
                    del tbl_old['samples']

                c = list(tbl_old.columns)
                c.remove(f"entity:{table_name}_id")
                c = [f"entity:{table_name}_id"] + c
                tbl_old = tbl_old[c]
                tbl_old = tbl_old.astype(str)

            return tbl_old, membership

        # Remove any `nan` values in a given dataframe.
        # `nan` values are caused by a parsing issue and are artifacts.
        def fix_nans(df, quiet=True):
            if not quiet: print("Replacing all `nan` values with empty strings: ")
            for c in df.columns.values:
                nan_types = ("nan", float('nan'))
                has_nan = False
                num_denaned = 0
                for n in nan_types:
                    if (sum(df[c] == n) > 0):
                        num_denaned += sum(df[c] == n)
                        df.loc[df[c] == n, c] = ""
                        has_nan = True
                if has_nan and not quiet:
                    print(f"\t{c}: {num_denaned}")

            if not quiet: print("Replacing numpy nan values...")
            if not quiet: print("Done.")
            return df.replace(np.nan, "")

        def _write_json_file_to_bucket(backup_folder_path, bucket, timestamp, namespace, workspace, json_object, file_base_name):
            # Write our dict to our bucket:
            blob = bucket.blob(f"{backup_folder_path}/{timestamp}_{namespace}_{workspace}_{file_base_name}.json.gz")
            with blob.open('wb') as f:
                f.write(gzip.compress(bytes(json.dumps(json_object, indent=4), 'utf-8')))


        ################################################################################

        # Get our entity types so we know what to dump:
        entity_types = fapi.list_entity_types(namespace, workspace).json()

        timestamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
        workspace_bucket = fapi.get_workspace(namespace, workspace).json()["workspace"]["bucketName"]
        backup_folder_path = f"backups/{timestamp}"

        # Create our timestamped backup bucket:
        storage_client = storage.Client()
        bucket = storage_client.bucket(workspace_bucket)

        ################################################################################

        # Iterate over entity types and dump each one to a separate TSV:
        print(f"Writing workspace entities to backup dir:")
        print(f"gs://{workspace_bucket}/{backup_folder_path}/tables")
        for et in entity_types:
            print(f"\t{et}")
            tbl, _ = load_table(namespace, workspace, et)
            tbl = fix_nans(tbl)
            table_name = f"{timestamp}_{namespace}_{workspace}_{et}.tsv.gz"

            # Write our table to our bucket:
            blob = bucket.blob(f"{backup_folder_path}/tables/{table_name}")

            with io.StringIO() as buf:
                tbl.to_csv(buf, sep="\t", index=False)
                with blob.open('wb') as f:
                    f.write(gzip.compress(bytes(buf.getvalue(), 'utf-8')))
        print('Done.')

        ################################################################################

        # Now backup the notebooks:
        print("Writing notebooks to backup dir:")
        print(f"gs://{workspace_bucket}/{backup_folder_path}/notebooks")
        for notebook_blob in storage_client.list_blobs(workspace_bucket, prefix='notebooks'):
            original_name = notebook_blob.name[notebook_blob.name.find("/")+1:]
            print(f"\t{original_name}")
            notebook_name = f"{timestamp}_{namespace}_{workspace}_{original_name}"
            blob = bucket.copy_blob(notebook_blob, bucket, new_name=f"{backup_folder_path}/notebooks/{notebook_name}")
        print("Done.")

        ################################################################################

        # Now store workspace attributes:
        ws_dict = lrma_workspace_utils._query_workspace(namespace, workspace, 2)
        ws_attributes_dict = ws_dict['attributes']

        # Write our dict to our bucket:
        print(f"Writing workspace attributes to backup dir:")
        print(f"gs://{workspace_bucket}/{backup_folder_path}")

        _write_json_file_to_bucket(f"{backup_folder_path}", bucket, timestamp, namespace, workspace, ws_attributes_dict,
                                   f"workspace_attributes")

        print(f"Done.")
        print()

        # and store the workspace metadata:
        del ws_dict['attributes']

        print(f"Writing workspace metadata to backup dir:")
        print(f"gs://{workspace_bucket}/{backup_folder_path}")

        # Write our dict to our bucket:
        _write_json_file_to_bucket(f"{backup_folder_path}", bucket, timestamp, namespace, workspace, ws_dict,
                                   f"workspace_metadata")

        print(f"Done.")

        ################################################################################

        # Now store workspace method (i.e. workflow) information:
        response = lrma_workspace_utils.retry_fiss_api_call('list_workspace_configs', 2, namespace, workspace, True)
        workflow_dict = response.json()

        # Sort our workflows in alphabetical order:
        workflow_dict = sorted(workflow_dict, key=lambda k: k['name'])

        print(f"Writing workflow high-level information to backup dir:")
        print(f"gs://{workspace_bucket}/{backup_folder_path}")

        # Write our dict to our bucket:
        _write_json_file_to_bucket(f"{backup_folder_path}", bucket, timestamp, namespace, workspace, workflow_dict,
                                       f"workflows")

        print(f'Done.')

        ################################################################################

        # Get the metadata, default inputs, and default outputs for each workflow:

        print(f"Writing workflow metadata, inputs, and outputs to backup dir:")
        print(f"gs://{workspace_bucket}/{backup_folder_path}/workflows")

        for workflow_name in [w["name"] for w in workflow_dict]:
            print(f"\t{workflow_name}")

            response = lrma_workspace_utils.retry_fiss_api_call('get_workspace_config', 2, namespace, workspace, namespace, workflow_name)
            metadata = response.json()

            # Write the workflow metadata:
            _write_json_file_to_bucket(f"{backup_folder_path}/workflows", bucket, timestamp, namespace, workspace, metadata,
                                       f"{workflow_name}_workflow_metadata")

            inputs = metadata['inputs']
            outputs = metadata['outputs']
            del metadata['inputs']
            del metadata['outputs']

            # Write the workflow inputs:
            _write_json_file_to_bucket(f"{backup_folder_path}/workflows", bucket, timestamp, namespace, workspace, inputs,
                                       f"{workflow_name}_workflow_inputs")
            # Write the workflow outputs:
            _write_json_file_to_bucket(f"{backup_folder_path}/workflows", bucket, timestamp, namespace, workspace, outputs,
                                       f"{workflow_name}_workflow_outputs")

        print(f"Done")

        ################################################################################

        # Finally store submissions:
        response = lrma_workspace_utils.retry_fiss_api_call('list_submissions', 2, namespace, workspace)
        submissions_dict = response.json()

        # Sort our submissions from most recent to least recent:
        submissions_dict = sorted(submissions_dict, key=lambda k: k['submissionDate'], reverse=True)

        print(f"Writing workspace job history to backup dir:")
        print(f"gs://{workspace_bucket}/{backup_folder_path}")

        _write_json_file_to_bucket(backup_folder_path, bucket, timestamp, namespace, workspace, submissions_dict, "workspace_job_history")

        print(f"Done.")

        ################################################################################

        # Track end time so we can calculate elapsed time for backup:
        END_TIME = time.time()

        ################################################################################

        import pytz
        now_utc = datetime.datetime.utcnow()
        timezone = pytz.timezone('America/New_York')
        now_et = now_utc.astimezone(timezone)
        time_string = now_et.strftime("%A %B %d at %H:%M:%S ET")
        print(f"Backup completed on {time_string}")
        print(f"Backup location: gs://{workspace_bucket}/{backup_folder_path}")

        # Write backup dir to a file so we can refer to it in our outputs:
        with open("~{backup_path_file}") as f:
            f.write(f"gs://{workspace_bucket}/{backup_folder_path}")

        print(f"Elapsed time: {END_TIME-START_TIME:2.2f}s")

        ################################################################################

        CODE
    >>>

    output {
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
