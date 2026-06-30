version 1.0

import "../../structs/Structs.wdl"

task UploadDataTable {

    meta {
        description: "Upload a Terra-formatted entity TSV into a workspace data table. Optionally require that every entity ID in the new data already exists in the destination table, failing (and listing the offenders on stderr) when any do not."
        volatile: true
    }

    parameter_meta {
        namespace:           "Terra billing project / namespace that owns the destination workspace."
        workspace:           "Name of the destination Terra workspace."
        entities_tsv:        "Terra-formatted entity TSV to upload. Its first column header must be `entity:<table_name>_id`."
        table_name:          "Name of the destination data table (entity type) to upload into."
        require_existing_id: "When true (default), pull the destination table first and verify every entity ID in `entities_tsv` already exists there; any IDs that do not are printed to stderr and the task fails before uploading. When false, upload unconditionally."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        String namespace
        String workspace
        File entities_tsv
        String table_name

        Boolean require_existing_id = true

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(size(entities_tsv, "GB"))

    command <<<
        set -euxo pipefail

        python3 <<CODE
        import sys

        import pandas as pd
        import firecloud.api as fapi

        namespace           = "~{namespace}"
        workspace           = "~{workspace}"
        table_name          = "~{table_name}"
        entities_tsv        = "~{entities_tsv}"
        require_existing_id = ~{true="True" false="False" require_existing_id}

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

        if require_existing_id:
            # Pull the destination table and collect the IDs that already exist.
            resp = fapi.get_entities(namespace, workspace, table_name)
            if resp.status_code != 200:
                sys.stderr.write(
                    f"ERROR: could not read destination table '{table_name}' from "
                    f"{namespace}/{workspace}: HTTP {resp.status_code}: {resp.text}\n"
                )
                sys.exit(1)
            existing_ids = set(e["name"] for e in resp.json())

            # Any new ID not already present is an error; list them all.
            missing = [i for i in new_ids if i not in existing_ids]
            if missing:
                sys.stderr.write(
                    f"ERROR: {len(missing)} of {len(new_ids)} sample ID(s) in {entities_tsv} are not "
                    f"present in destination table '{table_name}' ({namespace}/{workspace}):\n"
                )
                for m in missing:
                    sys.stderr.write(f"  {m}\n")
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
