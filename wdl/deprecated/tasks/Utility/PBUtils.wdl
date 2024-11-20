version 1.0

import "../../../structs/Structs.wdl"

task MakePerBarcodeDemultiplexingReports {

    meta {
        description: "Make per-barcode demultiplexing reports."
    }

    parameter_meta {
        report: "Input report."
        type: "Output file type (pdf or png)."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        File report
        String type = "png"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(report, "GB"))

    command <<<
        set -x

        grep '>' /Sequel_96_barcodes_v2.fasta | sed 's/>//' | while read -r line ; do
            Rscript /lima_report_detail.R ~{report} ~{type} $line

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
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
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
