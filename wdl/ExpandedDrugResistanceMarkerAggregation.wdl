version 1.0

import "tasks/Structs.wdl"
import "tasks/FunctionalAnnotation.wdl" as FUNK
import "tasks/Finalize.wdl" as FF

workflow ExpandedDrugResistanceMarkerExtraction {
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
    input {
        Array[String] sample_names
        Array[File] expanded_drug_res_markers

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + 10*ceil(size([sample_names, expanded_drug_res_markers], "GB"))

    command <<<
        set -euxo pipefail

        python3 <<CODE

        aggregated_report = "~{prefix}.expanded_drug_report_combined.tsv"

        sample_name_file = ~{write_lines(sample_names)}
        drug_res_file = ~{write_lines(expanded_drug_res_markers)}

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
                next(f)
                for line in f:
                    fields = line.strip().split()
                    marker = "_".join(fields[:3])
                    sample_marker_dict[sample_id][marker] = 1 if fields[4] == "hom" else -2

        # Make lists of the things we want to iterate over:
        marker_list_list = list([v.keys() for v in sample_marker_dict.values()])
        marker_list = sorted(list(set([m for ml in marker_list_list for m in ml])))

        sample_list = sorted([k for k in sample_marker_dict.keys()])

        # # Second pass to actually aggregate the markers:
        with open(aggregated_report, 'w') as f:
            # Write header:
            f.write("Marker\t" + "\t".join(sample_list))
            f.write("\n")
            for marker in marker_list:
                f.write(marker + "\t")
                for sample in sample_list:
                    f.write(str(sample_marker_dict[sample][marker]))
                    f.write("\t")
                f.write("\n")
        CODE

    >>>

    output {
        File combined_report = "~{prefix}.expanded_drug_report_combined.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/biocontainers/snpeff:5.1d--hdfd78af_0"
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
