version 1.0

import "../../structs/Structs.wdl"

task UploadDataTable {

    meta {
        description: "Upload a Terra-formatted entity TSV into a workspace data table. Optionally require that every entity ID in the new data already exists in the destination table, failing (and listing the offenders on stderr) when any do not."
        volatile: true
    }

    parameter_meta {
        namespace:           "Optional Terra billing project / namespace that owns the destination workspace. When omitted it is auto-detected (WORKSPACE_NAMESPACE env var, else the workspace whose Google project matches this VM)."
        workspace:           "Optional name of the destination Terra workspace. When omitted it is auto-detected (WORKSPACE_NAME env var, else the workspace whose Google project matches this VM)."
        entities_tsv:        "Terra-formatted entity TSV to upload. Its first column header must be `entity:<table_name>_id`."
        table_name:          "Name of the destination data table (entity type) to upload into."
        require_existing_id: "When true (default), pull the destination table first and verify every entity ID in `entities_tsv` already exists there; any IDs that do not are printed to stderr and the task fails before uploading. When false, upload unconditionally."
        columns_must_exist:  "When true (default), every column in `entities_tsv` must already exist in the destination table; any that do not are printed to stderr and the task fails before uploading. When false, new columns are allowed."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        String? namespace
        String? workspace
        File entities_tsv
        String table_name

        Boolean require_existing_id = true
        Boolean columns_must_exist = true

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(size(entities_tsv, "GB"))

    command <<<
        set -euxo pipefail

        python3 <<CODE
        import os
        import sys
        import urllib.request

        import pandas as pd
        import firecloud.api as fapi

        def gce_project():
            # Fetch this VM's Google project from the metadata server (no curl
            # dependency; the runtime image may not ship curl).
            req = urllib.request.Request(
                "http://metadata.google.internal/computeMetadata/v1/project/project-id",
                headers={"Metadata-Flavor": "Google"},
            )
            try:
                return urllib.request.urlopen(req, timeout=5).read().decode().strip()
            except Exception as e:
                sys.stderr.write("WARNING: could not read Google project from metadata server: {}\n".format(e))
                return ""

        table_name          = "~{table_name}"
        entities_tsv        = "~{entities_tsv}"
        require_existing_id = ~{true="True" false="False" require_existing_id}
        columns_must_exist  = ~{true="True" false="False" columns_must_exist}

        # namespace/workspace are optional; empty string means "auto-detect".
        namespace = "~{default='' namespace}"
        workspace = "~{default='' workspace}"

        def resolve_workspace(namespace, workspace):
            # 1. Explicit inputs win.
            if namespace and workspace:
                return namespace, workspace

            # 2. Terra-injected environment variables (set in interactive runtimes).
            namespace = namespace or os.environ.get("WORKSPACE_NAMESPACE", "")
            workspace = workspace or os.environ.get("WORKSPACE_NAME", "")
            if namespace and workspace:
                return namespace, workspace

            # 3. Derive from the Google project this VM runs in. In the current Terra
            #    model each workspace has its own Google project, so the project maps
            #    back to a unique workspace via the workspaces list.
            project = os.environ.get("GOOGLE_PROJECT", "") or gce_project()
            if not project:
                sys.stderr.write(
                    "ERROR: namespace/workspace not provided and the Google project could not be "
                    "determined (GOOGLE_PROJECT unset and metadata server unavailable); "
                    "cannot auto-detect the destination workspace.\n"
                )
                sys.exit(1)

            resp = fapi.list_workspaces(fields="workspace.namespace,workspace.name,workspace.googleProject")
            if resp.status_code != 200:
                sys.stderr.write(
                    f"ERROR: could not list workspaces to auto-detect the destination: "
                    f"HTTP {resp.status_code}: {resp.text}\n"
                )
                sys.exit(1)

            matches = sorted({
                (w["workspace"]["namespace"], w["workspace"]["name"])
                for w in resp.json()
                if w["workspace"].get("googleProject") == project
            })
            if len(matches) == 1:
                return matches[0]

            sys.stderr.write(
                f"ERROR: could not uniquely auto-detect the workspace for Google project "
                f"'{project}' (found {len(matches)} match(es)); pass namespace and workspace "
                f"explicitly.\n"
            )
            sys.exit(1)

        namespace, workspace = resolve_workspace(namespace, workspace)

        # Identify the entity-ID column. Terra entity TSVs name it 'entity:<table>_id'.
        id_col = f"entity:{table_name}_id"
        new_df = pd.read_csv(entities_tsv, sep="\t", dtype=str)
        if id_col not in new_df.columns:
            first = new_df.columns[0]
            if first.startswith("entity:") and first.endswith("_id"):
                id_col = first
            else:
                sys.stderr.write(
                    f"ERROR: expected an entity-ID column named '{id_col}' (or 'entity:<table>_id') "
                    f"in {entities_tsv}; found columns: {list(new_df.columns)}\n"
                )
                sys.exit(1)

        new_ids = [str(x) for x in new_df[id_col].tolist()]

        if require_existing_id or columns_must_exist:
            # Pull the destination table once to validate IDs and/or columns against it.
            resp = fapi.get_entities(namespace, workspace, table_name)
            if resp.status_code != 200:
                sys.stderr.write(
                    f"ERROR: could not read destination table '{table_name}' from "
                    f"{namespace}/{workspace}: HTTP {resp.status_code}: {resp.text}\n"
                )
                sys.exit(1)
            dest_entities = resp.json()

            fail = False

            if require_existing_id:
                existing_ids = set(e["name"] for e in dest_entities)
                # Any new ID not already present is an error; list them all.
                missing_ids = [i for i in new_ids if i not in existing_ids]
                if missing_ids:
                    sys.stderr.write(
                        f"ERROR: {len(missing_ids)} of {len(new_ids)} sample ID(s) in {entities_tsv} "
                        f"are not present in destination table '{table_name}' "
                        f"({namespace}/{workspace}):\n"
                    )
                    for m in missing_ids:
                        sys.stderr.write(f"  {m}\n")
                    fail = True

            if columns_must_exist:
                # Destination columns are the union of attribute keys across its rows.
                existing_cols = set()
                for e in dest_entities:
                    existing_cols.update(e.get("attributes", {}).keys())
                # The entity-ID column is not a data attribute, so exclude it.
                new_cols = [c for c in new_df.columns if c != id_col]
                missing_cols = [c for c in new_cols if c not in existing_cols]
                if missing_cols:
                    sys.stderr.write(
                        f"ERROR: {len(missing_cols)} column(s) in {entities_tsv} are not present in "
                        f"destination table '{table_name}' ({namespace}/{workspace}):\n"
                    )
                    for c in missing_cols:
                        sys.stderr.write(f"  {c}\n")
                    fail = True

            if fail:
                sys.exit(1)

        # All checks passed (or were skipped): upload the table.
        up = fapi.upload_entities_tsv(namespace, workspace, entities_tsv, model="flexible")
        if up.status_code not in (200, 201):
            sys.stderr.write(
                f"ERROR: upload to '{table_name}' ({namespace}/{workspace}) failed: "
                f"HTTP {up.status_code}: {up.text}\n"
            )
            sys.exit(1)

        print(f"Uploaded {len(new_ids)} row(s) to table '{table_name}' in {namespace}/{workspace}.")
        CODE
    >>>

    output {
        File upload_log = stdout()
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
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

task GetWorkspaceMetadata {

    meta {
        description: "Fetch the FULL Terra metadata for the workspace this task is running in (fapi.get_workspace). The workspace is auto-detected -- WORKSPACE_NAMESPACE/WORKSPACE_NAME env vars (set in Terra runtimes), else the workspace whose Google project matches this VM's project (metadata server) via fapi.list_workspaces -- unless namespace/workspace are passed. Emits the raw JSON, the workspace attributes as their own JSON, and the commonly-used scalar fields as typed outputs."
        volatile: true
    }

    parameter_meta {
        namespace:             "Optional Terra billing project / namespace. Auto-detected when omitted."
        workspace:             "Optional workspace name. Auto-detected when omitted."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        String? namespace
        String? workspace

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        python3 <<CODE
        import os, sys, json, urllib.request
        import firecloud.api as fapi

        def gce_project():
            # This VM's Google project from the metadata server (no curl dependency).
            req = urllib.request.Request(
                "http://metadata.google.internal/computeMetadata/v1/project/project-id",
                headers={"Metadata-Flavor": "Google"},
            )
            try:
                return urllib.request.urlopen(req, timeout=5).read().decode().strip()
            except Exception as e:
                sys.stderr.write("WARNING: could not read Google project from metadata server: {}\n".format(e))
                return ""

        namespace = "~{default='' namespace}"
        workspace = "~{default='' workspace}"

        def resolve_workspace(namespace, workspace):
            # 1. Explicit inputs win.
            if namespace and workspace:
                return namespace, workspace
            # 2. Terra-injected environment variables (interactive runtimes).
            namespace = namespace or os.environ.get("WORKSPACE_NAMESPACE", "")
            workspace = workspace or os.environ.get("WORKSPACE_NAME", "")
            if namespace and workspace:
                return namespace, workspace
            # 3. Derive from this VM's Google project (one project per workspace in Terra).
            project = os.environ.get("GOOGLE_PROJECT", "") or gce_project()
            if not project:
                sys.stderr.write("ERROR: no namespace/workspace and no Google project; cannot auto-detect workspace.\n")
                sys.exit(1)
            resp = fapi.list_workspaces(fields="workspace.namespace,workspace.name,workspace.googleProject")
            if resp.status_code != 200:
                sys.stderr.write(f"ERROR: list_workspaces HTTP {resp.status_code}: {resp.text}\n")
                sys.exit(1)
            matches = sorted({
                (w["workspace"]["namespace"], w["workspace"]["name"])
                for w in resp.json()
                if w["workspace"].get("googleProject") == project
            })
            if len(matches) == 1:
                return matches[0]
            sys.stderr.write(
                f"ERROR: could not uniquely auto-detect the workspace for Google project '{project}' "
                f"(found {len(matches)} match(es)); pass namespace and workspace explicitly.\n"
            )
            sys.exit(1)

        namespace, workspace = resolve_workspace(namespace, workspace)

        resp = fapi.get_workspace(namespace, workspace)
        if resp.status_code != 200:
            sys.stderr.write(f"ERROR: get_workspace {namespace}/{workspace} HTTP {resp.status_code}: {resp.text}\n")
            sys.exit(1)
        full  = resp.json()
        ws    = full.get("workspace", {}) or {}
        attrs = ws.get("attributes", {}) or {}

        # Full payload + the custom attributes as their own files.
        with open("workspace_metadata.json", "w") as f:
            json.dump(full, f, indent=2, sort_keys=True)
        with open("workspace_attributes.json", "w") as f:
            json.dump(attrs, f, indent=2, sort_keys=True)

        # Commonly-used scalars -> individual files for typed WDL outputs.
        def w(path, val):
            with open(path, "w") as f:
                f.write("" if val is None else str(val))

        w("namespace.txt",      ws.get("namespace", ""))
        w("name.txt",           ws.get("name", ""))
        w("workspace_id.txt",   ws.get("workspaceId", ""))
        w("bucket_name.txt",    ws.get("bucketName", ""))
        w("bucket_uri.txt",     ("gs://" + ws["bucketName"]) if ws.get("bucketName") else "")
        w("google_project.txt", ws.get("googleProject", ""))
        w("created_by.txt",     ws.get("createdBy", ""))
        w("created_date.txt",   ws.get("createdDate", ""))
        w("last_modified.txt",  ws.get("lastModified", ""))
        w("is_locked.txt",      "true" if ws.get("isLocked") else "false")
        w("access_level.txt",   full.get("accessLevel", ""))
        w("description.txt",    attrs.get("description", ""))
        ad = [g.get("membersGroupName", "") for g in (ws.get("authorizationDomain", []) or [])]
        with open("auth_domains.txt", "w") as f:
            f.write("\n".join([x for x in ad if x]))

        print(f"Fetched metadata for {namespace}/{workspace} (bucket gs://{ws.get('bucketName','?')}).")
        CODE
    >>>

    output {
        File          metadata_json         = "workspace_metadata.json"   # full fapi.get_workspace payload
        File          attributes_json       = "workspace_attributes.json" # just workspace.attributes (custom key-values)
        String        workspace_namespace   = read_string("namespace.txt")
        String        workspace_name        = read_string("name.txt")
        String        workspace_id          = read_string("workspace_id.txt")
        String        bucket_name           = read_string("bucket_name.txt")
        String        bucket                = read_string("bucket_uri.txt")  # gs://<bucket_name>
        String        google_project        = read_string("google_project.txt")
        String        created_by            = read_string("created_by.txt")
        String        created_date          = read_string("created_date.txt")
        String        last_modified         = read_string("last_modified.txt")
        Boolean       is_locked             = read_boolean("is_locked.txt")
        String        access_level          = read_string("access_level.txt")
        String        description           = read_string("description.txt")
        Array[String] authorization_domains = read_lines("auth_domains.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            10,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
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
