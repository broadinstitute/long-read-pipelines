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
        num_samples:         "Number of samples in the input VCF (optional; default: 100)."
        sample_names:        "Optional list of sample names to extract. If provided, every name must occur in input_vcf; the workflow fails when any is missing. Mutually exclusive with sample_name_list."
        sample_name_list:    "Optional file of sample names to extract, one per line. Mutually exclusive with sample_names."

        destination_table:   "Name of the destination Terra data table to write the per-sample rows into."
        joint_run_ID:        "Identifier of the joint-call run that produced input_vcf; recorded in the 97_Joint_Run_ID column of every row."

        namespace:           "Optional Terra namespace of the destination workspace (auto-detected when omitted)."
        workspace:           "Optional Terra workspace name (auto-detected when omitted)."
        require_existing_id: "When true, every sample ID must already exist in destination_table before upload. Default false (rows may be created)."
        columns_must_exist:  "When true, every column must already exist in destination_table before upload. Default false (columns may be created)."
    }

    input {
        File input_vcf
        File? input_vcf_index

        Int num_samples = 100

        Array[String]? sample_names
        File? sample_name_list

        String destination_table
        String joint_run_ID

        String? namespace
        String? workspace
        Boolean require_existing_id = false
        Boolean columns_must_exist = false
    }

    # 1. Split the joint VCF into one VCF (+ index) per sample.
    call Split.SplitMultiSampleVCFTask as SplitVCF {
        input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            num_samples_for_disk_size_scaling = num_samples,
            sample_names = sample_names,
            sample_name_list = sample_name_list
    }

    # 2. Build the Terra entity TSV, one row per per-sample VCF.
    call MakeEntitiesTsv {
        input:
            sample_vcfs = SplitVCF.output_vcfs,
            sample_vcf_indices = SplitVCF.output_vcf_indices,
            table_name = destination_table,
            joint_run_id = joint_run_ID
    }

    # 3. Upload the TSV into the destination data table.
    call Terra.UploadDataTable as Upload {
        input:
            namespace = namespace,
            workspace = workspace,
            entities_tsv = MakeEntitiesTsv.entities_tsv,
            table_name = destination_table,
            require_existing_id = require_existing_id,
            columns_must_exist = columns_must_exist
    }

    output {
        Array[File] sample_vcfs = SplitVCF.output_vcfs
        Array[File] sample_vcf_indices = SplitVCF.output_vcf_indices
        File entities_tsv = MakeEntitiesTsv.entities_tsv
        File upload_log = Upload.upload_log
    }
}

task MakeEntitiesTsv {

    meta {
        description: "Assemble a Terra entity TSV from a set of per-sample VCFs. Each row uses the sample name (parsed from the VCF file name) as the entity ID and records the sample's VCF, its index, and the joint-call run ID."
    }

    parameter_meta {
        sample_vcfs:        "Per-sample VCF files (one per sample); the sample name is parsed from each file's base name."
        sample_vcf_indices: "Per-sample VCF index files, matched to sample_vcfs by sample name."
        table_name:         "Destination data table name; used to form the `entity:<table_name>_id` header column."
        joint_run_id:       "Joint-call run ID recorded in the 97_Joint_Run_ID column of every row."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        Array[String] sample_vcfs
        Array[String] sample_vcf_indices
        String table_name
        String joint_run_id

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        cp "~{write_lines(sample_vcfs)}" vcfs.txt
        cp "~{write_lines(sample_vcf_indices)}" idxs.txt

        python3 <<CODE
        import os

        table_name   = "~{table_name}"
        joint_run_id = "~{joint_run_id}"

        def sample_of(path):
            b = os.path.basename(path)
            # Strip the longest matching VCF extension to recover the sample name.
            for ext in (".vcf.gz.tbi", ".vcf.gz", ".vcf.tbi", ".vcf", ".tbi"):
                if b.endswith(ext):
                    return b[:-len(ext)]
            return b

        with open("vcfs.txt") as fh:
            vcfs = [line.strip() for line in fh if line.strip()]
        with open("idxs.txt") as fh:
            idxs = [line.strip() for line in fh if line.strip()]

        vcf_by_sample = {sample_of(v): v for v in vcfs}
        idx_by_sample = {sample_of(i): i for i in idxs}

        with open("entities.tsv", "w") as out:
            out.write("entity:{}_id\tvcf\tvcf_index\t98_Run_ID\t97_Joint_Run_ID\n".format(table_name))
            for sample in sorted(vcf_by_sample):
                out.write("{s}\t{v}\t{i}\t{s}\t{j}\n".format(
                    s=sample,
                    v=vcf_by_sample[sample],
                    i=idx_by_sample.get(sample, ""),
                    j=joint_run_id,
                ))
        CODE
    >>>

    output {
        File entities_tsv = "entities.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            10,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
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
