version 1.0

##########################################################################################
## A workflow that processes P. falciparum SNP panels (read: barcodes) and calculates several
## metrics that are relevant to studying the epidemiology of the disease.
##
## This WDL calls a script written by Wes Wong and based on the following paper:
## https://doi.org/10.1093/pnasnexus/pgac187
##########################################################################################

import "tasks/Structs.wdl"
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF

workflow PanelProcessMalariaBarcodesForRh {
    input {

        # Unfortunately the easiest way to make this work would be to pass a spreadsheet into the script.
        # Because of how Terra is structured, this isn't really possible.
        # Instead, we pass each column and construct a spreadsheet in the task.
        # This way we can have one big table of the values in Terra and we don't have to make the data hard
        # to visualize.

        # High-level required info:
        String location_code
        File barcode_def_tsv

        # Spreadsheet data:
        Array[String] cc
        Array[String] ISO3
        Array[String] Year
        Array[String] Number_Text
        Array[String] Sample_Name
        Array[String] Raw_Name
        Array[String] Barcode_String
        Array[String] A1
        Array[String] B1
        Array[String] A2
        Array[String] B2
        Array[String] A3
        Array[String] B3
        Array[String] A4
        Array[String] B4
        Array[String] A5
        Array[String] B5
        Array[String] A6
        Array[String] B6
        Array[String] A7
        Array[String] B7
        Array[String] A8
        Array[String] B8
        Array[String] A9
        Array[String] B9
        Array[String] A10
        Array[String] B10
        Array[String] A11
        Array[String] B11
        Array[String] A12
        Array[String] B12
        Array[String] X
        Array[String] N
        Array[String] M_P
        Array[String] Delta_CT_Threshold
        Array[String] Adjusted_Het
        Array[String] mccoil_median

        String dir_prefix
        String gcs_out_root_dir

        Boolean DEBUG_MODE = false
    }

    parameter_meta {
        # TODO: FILL THIS IN
    }

    ####################################
    #     _____         _
    #    |_   _|_ _ ___| | _____
    #      | |/ _` / __| |/ / __|
    #      | | (_| \__ \   <\__ \
    #      |_|\__,_|___/_|\_\___/
    #
    ####################################

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_001_WdlExecutionStartTimestamp { input: }

    # Create an outdir:
    String outdir = if DEBUG_MODE then sub(gcs_out_root_dir, "/$", "") + "/PanelProcessMalariaBarcodesForRh/~{dir_prefix}/" + t_001_WdlExecutionStartTimestamp.timestamp_string else sub(gcs_out_root_dir, "/$", "") + "/PanelProcessMalariaBarcodesForRh/~{dir_prefix}"

    call ProcessBarcodeSpreadsheet {
        input:
            location_code = location_code,
            barcode_def_tsv = barcode_def_tsv,
            cc = cc,
            ISO3 = ISO3,
            Year = Year,
            Number_Text = Number_Text,
            Sample_Name = Sample_Name,
            Raw_Name = Raw_Name,
            Barcode_String = Barcode_String,
            A1 = A1,
            B1 = B1,
            A2 = A2,
            B2 = B2,
            A3 = A3,
            B3 = B3,
            A4 = A4,
            B4 = B4,
            A5 = A5,
            B5 = B5,
            A6 = A6,
            B6 = B6,
            A7 = A7,
            B7 = B7,
            A8 = A8,
            B8 = B8,
            A9 = A9,
            B9 = B9,
            A10 = A10,
            B10 = B10,
            A11 = A11,
            B11 = B11,
            A12 = A12,
            B12 = B12,
            X = X,
            N = N,
            M_P = M_P,
            Delta_CT_Threshold = Delta_CT_Threshold,
            Adjusted_Het = Adjusted_Het,
            mccoil_median = mccoil_median
    }
}

task ProcessBarcodeSpreadsheet {
    input {
        # High-level required info:
        String location_code
        File barcode_def_tsv

        # Spreadsheet data:
        Array[String] cc
        Array[String] ISO3
        Array[String] Year
        Array[String] Number_Text
        Array[String] Sample_Name
        Array[String] Raw_Name
        Array[String] Barcode_String
        Array[String] A1
        Array[String] B1
        Array[String] A2
        Array[String] B2
        Array[String] A3
        Array[String] B3
        Array[String] A4
        Array[String] B4
        Array[String] A5
        Array[String] B5
        Array[String] A6
        Array[String] B6
        Array[String] A7
        Array[String] B7
        Array[String] A8
        Array[String] B8
        Array[String] A9
        Array[String] B9
        Array[String] A10
        Array[String] B10
        Array[String] A11
        Array[String] B11
        Array[String] A12
        Array[String] B12
        Array[String] X
        Array[String] N
        Array[String] M_P
        Array[String] Delta_CT_Threshold
        Array[String] Adjusted_Het
        Array[String] mccoil_median

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    # Create a header for the inputs so we can generate a TSV:
    Array[String] header = ["cc", "ISO3", "Year", "Number_Text", "Sample_Name", "Raw_Name", "Barcode_String", "A1", "B1", "A2", "B2", "A3", "B3", "A4", "B4", "A5", "B5", "A6", "B6", "A7", "B7", "A8", "B8", "A9", "B9", "A10", "B10", "A11", "B11", "A12", "B12", "X", "N", "M_P", "Delta_CT_Threshold", "Adjusted_Het", "mccoil_median"]

    String input_tsv_path = "~{location_code}.tmp_input.tsv"

    command <<<
        set -euxo pipefail

        ## Generate the input TSV:
        tmp_input_tsv=~{write_tsv([header, cc, ISO3, Year, Number_Text, Sample_Name, Raw_Name, Barcode_String, A1, B1, A2, B2, A3, B3, A4, B4, A5, B5, A6, B6, A7, B7, A8, B8, A9, B9, A10, B10, A11, B11, A12, B12, X, N, M_P, Delta_CT_Threshold, Adjusted_Het, mccoil_median])}
        mv ${tmp_input_tsv} ~{input_tsv_path}

        ## Run the script:
        /python_scripts/process_barcode_data.py -b ~{barcode_def_tsv} -s ~{location_code} -f ~{input_tsv_path}
    >>>

    output {
        File summary_figure_svg = "~{location_code}_summary_figure.svg"
        File summary_figure_png = "~{location_code}_summary_figure.png"
        File summary_stats = "~{location_code}_summary.csv"
        File mono_barcode_stats = "~{location_code}_mono_barcodes.csv"
        File poly_barcode_stats = "~{location_code}_poly_barcodes.csv"

        File input_tsv = input_tsv_path
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-malaria:0.0.1"
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
