version 1.0

import "../../../structs/Structs.wdl"
import "./SplitMultiSampleVCF.wdl" as Split
import "../../../tasks/Utility/TerraUtils.wdl" as Terra

workflow SplitMultiSampleVCFToTerraTable {
    meta {
        description: "Split a joint-called multi-sample VCF into per-sample VCFs, assemble a Terra entity TSV (one row per sample) pointing at those VCFs, and upload it into a destination data table. By default the table and any missing columns are created on upload."
        author: "Jonn Smith"
    }

    parameter_meta {
        input_vcf:           "Joint-called multi-sample VCF file (can be compressed or uncompressed)."
        input_vcf_index:     "Index file for the input VCF (required if VCF is compressed)."
        disk_space_multiplier: "Multiplier applied to the input VCF size when sizing the split task's disk (optional; default: 4)."
        sample_names:        "Optional list of sample names to extract. If provided, every name must occur in input_vcf; the workflow fails when any is missing. Mutually exclusive with sample_name_list."
        sample_name_list:    "Optional file of sample names to extract, one per line. Mutually exclusive with sample_names."

        destination_table:   "Name of the destination Terra data table to write the per-sample rows into. Expected to follow '[TEST_]joint_results_<SET>_singlesample'. Also drives auto-derivation of joint_run_table/truth_table (unless those are given): joint_run_table is this name without the trailing '_singlesample' (TEST_ preserved), and truth_table is 'truth_<SET>'."
        joint_run_ID:        "Entity name of the joint-call run that produced input_vcf; every row's 97_Joint_Run_ID column links to it."
        joint_run_table:     "Optional. Entity type (data table) that joint_run_ID lives in. When omitted, derived as 'TEST_joint_results_<SET>'. Existence of this table (and the joint_run_ID entity in it) is verified before upload."
        truth_table:         "Optional. Truth-set entity type (data table). When omitted, derived as 'truth_<SET>'. Existence of this table (and each linked truth entity) is verified before upload. Each row's Truth column links to the truth entity named by the sample prefix (row ID up to the first underscore)."

        namespace:           "Optional Terra namespace of the destination workspace (auto-detected when omitted)."
        workspace:           "Optional Terra workspace name (auto-detected when omitted)."
        require_existing_id: "When true, every sample ID must already exist in destination_table before upload. Default false (rows may be created)."
        columns_must_exist:  "When true, every column must already exist in destination_table before upload. Default false (columns may be created)."
    }

    input {
        File input_vcf
        File? input_vcf_index

        Int disk_space_multiplier = 4

        Array[String]? sample_names
        File? sample_name_list

        String destination_table
        String joint_run_ID
        String? joint_run_table
        String? truth_table

        String? namespace
        String? workspace
        Boolean require_existing_id = false
        Boolean columns_must_exist = false
    }

    # Auto-derive the reference-target table names from the destination table,
    # unless the caller supplied them explicitly. Both are existence-checked below.
    # Convention: destination is "[TEST_]joint_results_<SET>_singlesample". The
    # joint-run table is the same name without the trailing "_singlesample" (the
    # TEST_ prefix, if any, is preserved). The truth set token is that with the
    # optional TEST_ and joint_results_ prefix stripped; truth tables are never
    # TEST_-prefixed. E.g. TEST_joint_results_pfcrosses_v2_singlesample ->
    # joint-run TEST_joint_results_pfcrosses_v2, truth truth_pfcrosses_v2.
    String joint_run_table_derived = sub(destination_table, "_singlesample$", "")
    String truth_set = sub(joint_run_table_derived, "^(TEST_)?joint_results_", "")
    String resolved_joint_run_table = select_first([joint_run_table, joint_run_table_derived])
    String resolved_truth_table = select_first([truth_table, "truth_" + truth_set])

    # 1. Split the joint VCF into one VCF (+ index) per sample.
    call Split.SplitMultiSampleVCFTask as t_01_SplitVCF {
        input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            disk_space_multiplier = disk_space_multiplier,
            sample_names = sample_names,
            sample_name_list = sample_name_list
    }

    # The split task delocalizes each VCF and its index together; split the bundle
    # back into typed arrays (co-location preserved so each .tbi stays beside its VCF).
    scatter (f in t_01_SplitVCF.output_files) {
        Boolean is_index = (basename(f, ".tbi") != basename(f))
        if (!is_index) { File split_vcf = f }
        if (is_index)  { File split_idx = f }
    }
    Array[File] colocated_vcfs = select_all(split_vcf)
    Array[File] colocated_vcf_indices = select_all(split_idx)

    # 2. Build the Terra entity TSV, one row per per-sample VCF.
    call MakeEntitiesTsv as t_02_MakeEntitiesTsv {
        input:
            sample_vcfs = colocated_vcfs,
            sample_vcf_indices = colocated_vcf_indices,
            table_name = destination_table,
            joint_run_id = joint_run_ID,
            joint_run_table = resolved_joint_run_table,
            truth_table = resolved_truth_table
    }

    # 3. Sanity check: every link target must exist in Terra. Pull the joint-run
    #    and truth tables and confirm the joint_run_ID entity and each linked truth
    #    entity are present; fail (listing all offenders) before any upload. On
    #    success the TSV is passed through, so Upload is gated on this check.
    call ValidateTerraLinks as t_03_ValidateTerraLinks {
        input:
            namespace = namespace,
            workspace = workspace,
            entities_tsv = t_02_MakeEntitiesTsv.entities_tsv,
            sample_ids = t_02_MakeEntitiesTsv.sample_ids,
            joint_run_table = resolved_joint_run_table,
            joint_run_id = joint_run_ID,
            truth_table = resolved_truth_table
    }

    # 4. Upload the validated TSV into the destination data table.
    call Terra.UploadDataTable as t_04_Upload {
        input:
            namespace = namespace,
            workspace = workspace,
            entities_tsv = t_03_ValidateTerraLinks.validated_tsv,
            table_name = destination_table,
            require_existing_id = require_existing_id,
            columns_must_exist = columns_must_exist
    }

    output {
        Array[File] sample_vcfs = colocated_vcfs
        Array[File] sample_vcf_indices = colocated_vcf_indices
        File entities_tsv = t_02_MakeEntitiesTsv.entities_tsv
        File validation_log = t_03_ValidateTerraLinks.validation_log
        File upload_log = t_04_Upload.upload_log
    }
}

task MakeEntitiesTsv {

    meta {
        description: "Assemble a Terra entity TSV from a set of per-sample VCFs. Each row uses the sample name (parsed from the VCF file name) as the entity ID and records the sample's VCF and index, the run ID (98_Run_ID, same value as the entity ID), an entity reference (link) to the joint-call run, and an entity reference (link) to the matching sample in the truth set."
    }

    parameter_meta {
        sample_vcfs:        "Per-sample VCF files (one per sample); the sample name is parsed from each file's base name."
        sample_vcf_indices: "Per-sample VCF index files, matched to sample_vcfs by sample name."
        table_name:         "Destination data table name; used to form the `entity:<table_name>_id` header column."
        joint_run_id:       "Joint-call run entity name; every row's 97_Joint_Run_ID column links to it."
        joint_run_table:    "Entity type (data table) that joint_run_id lives in; the 97_Joint_Run_ID reference targets this table."
        truth_table:        "Truth-set entity type (data table). Each row's Truth column links to the truth entity whose name is the sample prefix (the row ID up to the first underscore)."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        Array[String] sample_vcfs
        Array[String] sample_vcf_indices
        String table_name
        String joint_run_id
        String joint_run_table
        String truth_table

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        cp "~{write_lines(sample_vcfs)}" vcfs.txt
        cp "~{write_lines(sample_vcf_indices)}" idxs.txt

        python3 <<CODE
        import os
        import json

        table_name      = "~{table_name}"
        joint_run_id    = "~{joint_run_id}"
        joint_run_table = "~{joint_run_table}"
        truth_table     = "~{truth_table}"

        def sample_of(path):
            b = os.path.basename(path)
            # Strip the longest matching VCF extension to recover the sample name.
            for ext in (".vcf.gz.tbi", ".vcf.gz", ".vcf.tbi", ".vcf", ".tbi"):
                if b.endswith(ext):
                    return b[:-len(ext)]
            return b

        def ref(entity_type, entity_name):
            # Terra entity reference; a JSON object cell is imported as a link
            # (clickable reference) by the flexible TSV importer.
            return json.dumps({"entityType": entity_type, "entityName": entity_name},
                              separators=(",", ":"))

        with open("vcfs.txt") as fh:
            vcfs = [line.strip() for line in fh if line.strip()]
        with open("idxs.txt") as fh:
            idxs = [line.strip() for line in fh if line.strip()]

        vcf_by_sample = {sample_of(v): v for v in vcfs}
        idx_by_sample = {sample_of(i): i for i in idxs}

        # Every row links to the same joint-call run entity.
        joint_ref = ref(joint_run_table, joint_run_id)

        samples = sorted(vcf_by_sample)

        with open("entities.tsv", "w") as out:
            out.write("entity:{}_id\tvcf\tvcf_index\t98_Run_ID\t97_Joint_Run_ID\tTruth\n".format(table_name))
            for sample in samples:
                # Truth entity name is the sample prefix: the row ID up to the first underscore.
                truth_name = sample.split("_", 1)[0]
                out.write("{s}\t{v}\t{i}\t{s}\t{j}\t{t}\n".format(
                    s=sample,
                    v=vcf_by_sample[sample],
                    i=idx_by_sample.get(sample, ""),
                    j=joint_ref,
                    t=ref(truth_table, truth_name),
                ))

        # Emit the row IDs so downstream validation reuses the exact same sample set.
        with open("sample_ids.txt", "w") as out:
            for sample in samples:
                out.write(sample + "\n")
        CODE
    >>>

    output {
        File entities_tsv = "entities.tsv"
        Array[String] sample_ids = read_lines("sample_ids.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            10,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "python:3.11-slim"
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

task ValidateTerraLinks {

    meta {
        description: "Verify that every entity reference the TSV will write actually exists in Terra: the joint-run table must contain the joint_run_id entity, and the truth table must contain the truth entity (sample prefix) for every sample. Fails (listing all offenders on stderr) before any upload. On success the input TSV is passed through unchanged so the uploader is gated on this check."
        volatile: true
    }

    parameter_meta {
        namespace:        "Optional Terra namespace of the workspace holding the link-target tables (auto-detected when omitted)."
        workspace:        "Optional Terra workspace name (auto-detected when omitted)."
        entities_tsv:     "The entity TSV to be uploaded; returned unchanged as validated_tsv when all checks pass."
        sample_ids:       "Destination row IDs; the truth entity checked for each is its sample prefix (ID up to the first underscore)."
        joint_run_table:  "Entity type that must contain joint_run_id (the 97_Joint_Run_ID link target)."
        joint_run_id:     "Entity name that must exist in joint_run_table."
        truth_table:      "Entity type that must contain each sample's truth entity (the Truth link target)."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        String? namespace
        String? workspace
        File entities_tsv
        Array[String] sample_ids
        String joint_run_table
        String joint_run_id
        String truth_table

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        cp "~{write_lines(sample_ids)}" sample_ids.txt

        python3 <<CODE
        import os
        import sys
        import urllib.request

        import firecloud.api as fapi

        joint_run_table = "~{joint_run_table}"
        joint_run_id    = "~{joint_run_id}"
        truth_table     = "~{truth_table}"

        namespace = "~{default='' namespace}"
        workspace = "~{default='' workspace}"

        def gce_project():
            req = urllib.request.Request(
                "http://metadata.google.internal/computeMetadata/v1/project/project-id",
                headers={"Metadata-Flavor": "Google"},
            )
            try:
                return urllib.request.urlopen(req, timeout=5).read().decode().strip()
            except Exception as e:
                sys.stderr.write("WARNING: could not read Google project from metadata server: {}\n".format(e))
                return ""

        def resolve_workspace(namespace, workspace):
            if namespace and workspace:
                return namespace, workspace
            namespace = namespace or os.environ.get("WORKSPACE_NAMESPACE", "")
            workspace = workspace or os.environ.get("WORKSPACE_NAME", "")
            if namespace and workspace:
                return namespace, workspace
            project = os.environ.get("GOOGLE_PROJECT", "") or gce_project()
            if not project:
                sys.stderr.write(
                    "ERROR: namespace/workspace not provided and the Google project could not be "
                    "determined; cannot locate the link-target tables.\n"
                )
                sys.exit(1)
            resp = fapi.list_workspaces(fields="workspace.namespace,workspace.name,workspace.googleProject")
            if resp.status_code != 200:
                sys.stderr.write("ERROR: could not list workspaces: HTTP {}: {}\n".format(resp.status_code, resp.text))
                sys.exit(1)
            matches = sorted({
                (w["workspace"]["namespace"], w["workspace"]["name"])
                for w in resp.json()
                if w["workspace"].get("googleProject") == project
            })
            if len(matches) == 1:
                return matches[0]
            sys.stderr.write(
                "ERROR: could not uniquely auto-detect the workspace for Google project '{}' "
                "(found {} match(es)); pass namespace and workspace explicitly.\n".format(project, len(matches))
            )
            sys.exit(1)

        def entity_names(namespace, workspace, table):
            # Returns the set of entity names in 'table', or None if the table is unreadable.
            resp = fapi.get_entities(namespace, workspace, table)
            if resp.status_code != 200:
                sys.stderr.write(
                    "ERROR: could not read table '{}' from {}/{}: HTTP {}: {}\n".format(
                        table, namespace, workspace, resp.status_code, resp.text)
                )
                return None
            return set(e["name"] for e in resp.json())

        namespace, workspace = resolve_workspace(namespace, workspace)

        with open("sample_ids.txt") as fh:
            sample_ids = [line.strip() for line in fh if line.strip()]

        fail = False

        # 1. Joint-run link target: table exists and contains joint_run_id.
        joint_names = entity_names(namespace, workspace, joint_run_table)
        if joint_names is None:
            fail = True
        elif joint_run_id not in joint_names:
            sys.stderr.write(
                "ERROR: joint-run entity '{}' not found in table '{}' ({}/{}).\n".format(
                    joint_run_id, joint_run_table, namespace, workspace)
            )
            fail = True

        # 2. Truth link targets: table exists and contains every sample's truth entity.
        truth_names = entity_names(namespace, workspace, truth_table)
        if truth_names is None:
            fail = True
        else:
            wanted = sorted({s.split("_", 1)[0] for s in sample_ids})
            missing = [t for t in wanted if t not in truth_names]
            if missing:
                sys.stderr.write(
                    "ERROR: {} truth entit(y/ies) missing from table '{}' ({}/{}):\n".format(
                        len(missing), truth_table, namespace, workspace)
                )
                for t in missing:
                    sys.stderr.write("  {}\n".format(t))
                fail = True

        if fail:
            sys.exit(1)

        print("All link targets exist: joint-run '{}' in '{}', and {} truth entit(y/ies) in '{}'.".format(
            joint_run_id, joint_run_table, len({s.split('_', 1)[0] for s in sample_ids}), truth_table))
        CODE

        # Reached only when every check passed; pass the TSV through to gate the upload.
        cp "~{entities_tsv}" validated.tsv
    >>>

    output {
        File validated_tsv = "validated.tsv"
        File validation_log = stdout()
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            10,
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
