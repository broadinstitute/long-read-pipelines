version 1.0

import "../../structs/Structs.wdl"

task CCSLima {
    meta {
        desciption: "A task for demultiplexing/un-barcoding on-instrument CCSed SMRTCell"
    }
    input {
        File bam
        File barcode_file
        String movie
        String barcode_design

        Array[String] custom_options

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        barcode_design: "Please see https://lima.how/barcode-design.html"
    }

    Int disk_size = 4*ceil(size(bam, "GB")) + 1

    command <<<
        set -euxo pipefail

        # see https://lima.how/get-started.html#example-executions
        ccs_options=("--ccs")
        if [[ ~{bam} == *hifi_reads.bam ]]; then
            if echo ~{barcode_design} | grep -qi "asymmetr"; then
                barcode_design_arg='ASYMMETRIC'
            elif echo ~{barcode_design} | grep -qiE "(symmetr|tail)"; then
                barcode_design_arg='SYMMETRIC'
            else
                echo "Unrecognized barcode design for hifi reads: ~{barcode_design}." && exit 1
            fi
            ccs_options=("--hifi-preset" "${barcode_design_arg}")
        elif echo ~{barcode_design} | grep -qi "asymmetric"; then
            ccs_options=("--ccs" "--different")
        elif echo ~{barcode_design} | grep -qiE "(symmetric|tail)"; then
            ccs_options=("--ccs" "--same")
        elif echo ~{barcode_design} | grep -qi "single" ; then
            ccs_options=("--ccs" "--min-score" "80" "--single-side")
        else
            echo "Unrecognized barcode design: ~{barcode_design}" && exit 1
        fi

        custom_options=(~{sep=' ' custom_options})

        # see https://lima.how/faq/how-to-run.html#how-to-run for peek-guess
        lima \
            "${custom_options[@]}" \
            "${ccs_options[@]}" \
            ~{bam} \
            ~{barcode_file} \
            ~{movie}.bam

        find . -type f -exec ls -lah {} \; | sort -k5,5 -h

        mkdir -p lima_reports
        mv ~{movie}.lima* lima_reports
        mv ~{movie}.json ~{movie}.consensusreadset.xml lima_reports

        mkdir -p unbarcoded
        mv ~{movie}.unbarcoded.* unbarcoded

        ls ~{movie}.*.bam | grep -F 'bc' > bams.list
        awk -F '.' '{print $2}' bams.list > barcodes.list
        cat bams.list
        cat barcodes.list
        if echo ~{barcode_design} | grep -qiE "(single|asymmetr)"; then
            echo
        elif echo ~{barcode_design} | grep -qiE "(symmetr|tail)"; then
            mv barcodes.list tmp.txt
            awk -F '-' '{print $1}' tmp.txt > barcodes.list
            cat barcodes.list
        else
            echo "are you kidding me?" && exit 1
        fi
        paste barcodes.list bams.list > bc_2_bam.tsv
    >>>

    output {
        Array[String] barcodes_list = read_lines("barcodes.list")

        Map[String, String] barcode_2_bam = read_map("bc_2_bam.tsv")

        Array[File] unbarcoded_files = glob("unbarcoded/*")

        Array[File] demux_bams = glob("~{movie}.*.bam")
        Array[File] demux_pbis = glob("~{movie}.*.pbi")
        Array[File] demux_xmls = glob("~{movie}.*.xml")

        File lima_counts = "lima_reports/~{movie}.lima.counts"
        File? lima_guess = "lima_reports/~{movie}.lima.guess"
        File lima_report = "lima_reports/~{movie}.lima.report"
        File lima_summary = "lima_reports/~{movie}.lima.summary"
        File lima_json = "lima_reports/~{movie}.json"
        File lima_xml = "lima_reports/~{movie}.consensusreadset.xml"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.36"
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

task ParseLimaOptionsFile {
    meta {
        desciption: "Given options file, parse it into an array that can be fed to lima directly"
    }
    input {
        File lima_options_file
    }
    parameter_meta {
        lima_options_file: "A 2-col TSV file with lima option long name in the 1st col, and value in the 2nd col. For flag-like options, the 2nd column should be set to 'true'. Omit if default values are desired."
    }

    command <<<

        awk '{print "--"$1}' ~{lima_options_file} > names.txt
        awk '{print $2}' ~{lima_options_file} > values.txt

        paste -d '\n' names.txt values.txt | grep -v "^true$"
    >>>

    output {
        Array[String] parsed_options = read_lines(stdout())
    }

    runtime {
        disks: "local-disk 50 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task LocateBarcodeSpecificFolders {
    meta {
        desciption: "Generates a 2-col TSV where the 1st col is the barcode name, and 2nd is the GCS \'folder\' holding the de-barcoded files."
    }
    input {
        Array[String] barcodes_list
        Array[String] finalized_dir_for_each_barcode
    }

    command <<<
        set -eux
        paste <(cat ~{write_lines(barcodes_list)}) \
              <(cat ~{write_lines(finalized_dir_for_each_barcode)}) \
              > barcode_2_folder.tsv
    >>>

    output {
        File barcode_2_dir = "barcode_2_folder.tsv"
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

###########################################
workflow GatherLimaAndCustomDemultiplexingMetrics {
    meta {
        desciption: "Generate custom performance reports by parsing lima.report, and gather them together with the raw metric files produced by lima itself."
    }
    input {
        File lima_report
        File lima_counts
        File? lima_guess
        File lima_summary
        File lima_json
        File lima_xml

        File barcodes_fasta
    }

    # make reports on demultiplexing
    call MakeSummarizedCustomDemultiplexingReport as SummarizedDemuxReportPNG { input: lima_report =lima_report }
    call MakeDetailedCustomDemultiplexingReport as DetailedDemuxReportPNG { input: lima_report =lima_report, type="png" }
    call MakeDetailedCustomDemultiplexingReport as DetailedDemuxReportPDF { input: lima_report =lima_report, type="pdf" }
    call MakePerBarcodeCustomDemultiplexingReports as PerBarcodeDetailedDemuxReportPNG { input: lima_report =lima_report, type="png", barcodes_fasta = barcodes_fasta }
    call MakePerBarcodeCustomDemultiplexingReports as PerBarcodeDetailedDemuxReportPDF { input: lima_report =lima_report, type="pdf", barcodes_fasta = barcodes_fasta }

    # gather all these reports, PDFs, PNGs, create zip, then finalize
    call GatherLimaReports {
        input:
            lima_raw_metric_files = select_all([lima_report, lima_counts, lima_guess, lima_summary, lima_json, lima_xml]),
            summary_figures = SummarizedDemuxReportPNG.report_files,
            detailed_pngs = DetailedDemuxReportPNG.report_files,
            detailed_pdfs = DetailedDemuxReportPDF.report_files,
            per_barcode_pngs = PerBarcodeDetailedDemuxReportPNG.report_files,
            per_barcode_pdfs = PerBarcodeDetailedDemuxReportPDF.report_files
    }

    output {
        File all_metrics_gz = GatherLimaReports.tar_gz_all_lima_metrics
    }
}

task GatherLimaReports {
    meta {
        desciption: "A trivial task for grouping together the raw and custom lima reports."
    }
    input {
        Array[File] lima_raw_metric_files
        Array[File] summary_figures
        Array[File] detailed_pngs
        Array[File] detailed_pdfs
        Array[File] per_barcode_pngs
        Array[File] per_barcode_pdfs
    }

    command <<<

        mkdir -p lima_metrics

        cp -r ~{sep=' ' lima_raw_metric_files}  lima_metrics/raw_metric_files/

        cp -r ~{sep=' ' summary_figures}  lima_metrics/summary_figures/
        cp -r ~{sep=' ' detailed_pngs}  lima_metrics/detailed_pngs/
        cp -r ~{sep=' ' detailed_pdfs}  lima_metrics/detailed_pdfs/
        cp -r ~{sep=' ' per_barcode_pngs}  lima_metrics/per_barcode_pngs/
        cp -r ~{sep=' ' per_barcode_pdfs}  lima_metrics/per_barcode_pdfs/

        tar -zcvf lima_metrics.tar.gz -C lima_metrics/ .
    >>>

    output {
        File tar_gz_all_lima_metrics = "lima_metrics.tar.gz"
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

###########################################
task MakeSummarizedCustomDemultiplexingReport {
    meta {
        desciption: "Custom report on demultiplexing; summarized."
    }
    input {
        File lima_report

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 100

    command <<<
        set -euxo pipefail

        Rscript /lima_report_summary.R ~{lima_report}
    >>>

    output {
        Array[File] report_files = glob("summary_*.png")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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

task MakeDetailedCustomDemultiplexingReport {
    meta {
        desciption: "Custom report on demultiplexing; in detail."
    }
    input {
        File lima_report
        String type = "png"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 100

    command <<<
        set -euxo pipefail

        Rscript /lima_report_detail.R ~{lima_report} ~{type}
    >>>

    output {
        Array[File] report_files = glob("detail_*~{type}")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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

task MakePerBarcodeCustomDemultiplexingReports {
    meta {
        desciption: "Custom report on demultiplexing; per-barcode detected."
    }
    input {
        File lima_report
        String type = "png"

        File barcodes_fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 100

    command <<<
        set -x

        grep '>' ~{barcodes_fasta} | sed 's/>//' | while read -r line ; do
            Rscript /lima_report_detail.R ~{lima_report} ~{type} $line

            if [ -f "detail_hq_length_hist_barcoded_or_not.~{type}" ]; then
                for f in detail_*; do mv $f $line.$f; done
            fi
        done
    >>>

    output {
        Array[File] report_files = glob("*.~{type}")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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
