version 1.0

import "../../../structs/Structs.wdl"
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ExpandedDrugResistanceMarkerAggregation {

    meta {
        author: "Jonn Smith"
        description: "Combine multiple expanded drug resistance marker reports together for a set of samples."
    }

    parameter_meta {
        sample_names: "Array of sample names corresponding to each file in `expanded_drug_res_markers`."
        expanded_drug_res_markers: "Array of per-sample drug resistance marker reports to aggregate."
        out_file_prefix: "Prefix to use for output files."
        dir_prefix: "Directory prefix to use for finalized location."
        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs."
    }

    input {
        Array[String] sample_names
        Array[File] expanded_drug_res_markers

        String out_file_prefix

        String dir_prefix
        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ExpandedDrugResistanceMarkerAggregation/~{dir_prefix}"

    call CombineExpandedDrugResistanceMarkers {
        input:
            expanded_drug_res_markers = expanded_drug_res_markers,
            sample_names = sample_names,
            prefix = out_file_prefix,
    }


    # Finalize data
    String dir = outdir + "/reports"

    call FF.FinalizeToFile as FinalizeDRReportAllMarkers { input: outdir = dir, file = CombineExpandedDrugResistanceMarkers.combined_report }

    output {
        File combined_expanded_markers = FinalizeDRReportAllMarkers.gcs_path
    }
}

task CombineExpandedDrugResistanceMarkers {

    meta {
        author: "Jonn Smith"
        description: "Combine multiple expanded drug resistance marker reports together for a set of samples."
    }

    parameter_meta {
        sample_names: "Array of sample names corresponding to each file in `expanded_drug_res_markers`."
        expanded_drug_res_markers: {
            description: "Array of per-sample drug resistance marker reports to aggregate.",
            localization_optional: true
        }
        prefix: "Prefix to use for output files."
        runtime_attr_override: "Override for default runtime attributes."
    }

    input {
        Array[String] sample_names
        Array[File] expanded_drug_res_markers

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + 10*ceil(size(expanded_drug_res_markers, "GB"))

    command <<<
        set -euxo pipefail

        # Copy the files from the cloud to this machine here so it goes faster:
        mkdir -p expanded_drug_reports
        cd expanded_drug_reports
        remote_sample_files=~{write_lines(expanded_drug_res_markers)}
        cat ${remote_sample_files} | gsutil -m cp -I .

        # Create local file list:
        cat ${remote_sample_files} | sed -e 's@^.*/@@' -e "s@^@$(pwd)/@g" > local_sample_files.txt

        cd ..

        python3 <<CODE
        import gzip
        from collections import defaultdict

        aggregated_report = "~{prefix}.expanded_drug_report_combined.tsv.gz"

        sample_name_file = "~{write_lines(sample_names)}"
        drug_res_file = "expanded_drug_reports/local_sample_files.txt"

        # Read in the sample names:
        sample_names = []
        with open(sample_name_file, 'r') as f:
            for line in f:
                sample_names.append(line.strip())

        # Read in the drug resistance files:
        drug_res_files = []
        with open(drug_res_file, 'r') as f:
            for line in f:
                drug_res_files.append(line.strip())

        # Make sure the sample name and drug resistance files are the same length:
        assert len(sample_names) == len(drug_res_files), "Sample names and drug resistance files are not the same length!"

        # First pass to get the list of markers:
        markers = set()
        sample_marker_dict = defaultdict(lambda: defaultdict(int))
        for i, expanded_drug_report in enumerate(drug_res_files):
            sample_id = sample_names[i]
            with open(expanded_drug_report, 'r') as f:
                # NOTE: Header should be the same for all samples, so we can reuse it later
                header = next(f).strip().split("\t")

                gt_field = header.index("GT")
                error_field = header.index("ERRORS / WARNINGS / INFO")

                for line in f:
                    fields = [fld.strip() for fld in line.split("\t")]
                    gt = fields[gt_field]

                    marker = tuple([item for i, item in enumerate(fields) if i not in (gt_field, error_field)])

                    sample_marker_dict[sample_id][marker] = 1 if (gt=="1/1" or gt=="1|1") else -2

        # Make lists of the things we want to iterate over:
        marker_list_list = list([v.keys() for v in sample_marker_dict.values()])
        marker_list = sorted(list(set([m for ml in marker_list_list for m in ml])))

        sample_list = sorted([k for k in sample_marker_dict.keys()])

        # # Second pass to actually aggregate the markers:
        out_header = [item for i, item in enumerate(header) if i not in (gt_field, error_field)]
        with gzip.open(aggregated_report, 'wt') as f:
            # Write header:
            f.write("\t".join(out_header) + "\t" + "\t".join(sample_list))
            f.write("\n")
            for marker in marker_list:
                # Markers are now tuples, so we need to write each element first:
                f.write("\t".join(marker))
                f.write("\t")
                for i, sample in enumerate(sample_list):
                    f.write(str(sample_marker_dict[sample][marker]))
                    if i != (len(sample_list)-1):
                        f.write("\t")
                f.write("\n")
        CODE

    >>>

    output {
        File combined_report = "~{prefix}.expanded_drug_report_combined.tsv.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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
