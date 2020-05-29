version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.14/wdl/tasks/Structs.wdl"

task NanoPlotFromSummary {
    input {
        Array[File] summary_files
        String plot_type

        RuntimeAttr? runtime_attr_override
    }

    Int num_cpus = 4
    Int disk_size = 2*ceil(size(summary_files, "GB"))

    command <<<
        set -euxo pipefail

        /usr/bin/time -v python3 /usr/local/lib/python3.7/site-packages/nanoplot/NanoPlot.py -t ~{num_cpus} -c orangered --N50 -f ~{plot_type} --summary ~{sep=' ' summary_files}
    >>>

    output {
        Array[File] plots = glob("*.~{plot_type}")
        File stats = "NanoStats.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          "~{num_cpus}",
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/biocontainers/nanoplot:1.28.0--py_0"
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

task NanoPlotFromAlignedBam {
    input {
        File unaligned_bam
        String plot_type

        RuntimeAttr? runtime_attr_override
    }

    Int num_cpus = 4
    Int disk_size = 2*ceil(size(unaligned_bam, "GB"))

    command <<<
        set -euxo pipefail

        python3 /usr/local/lib/python3.7/site-packages/nanoplot/NanoPlot.py -t ~{num_cpus} -c orangered -f ~{plot_type} --ubam ~{unaligned_bam}
    >>>

    output {
        File ActivePores_Over_Time = "ActivePores_Over_Time.~{plot_type}"
        File ActivityMap_ReadsPerChannel = "ActivityMap_ReadsPerChannel.~{plot_type}"
        File CumulativeYieldPlot_Gigabases = "CumulativeYieldPlot_Gigabases.~{plot_type}"
        File CumulativeYieldPlot_NumberOfReads = "CumulativeYieldPlot_NumberOfReads.~{plot_type}"
        File HistogramReadlength = "HistogramReadlength.~{plot_type}"
        File LengthvsQualityScatterPlot_dot = "LengthvsQualityScatterPlot_dot.~{plot_type}"
        File LengthvsQualityScatterPlot_kde = "LengthvsQualityScatterPlot_kde.~{plot_type}"
        File LogTransformed_HistogramReadlength = "LogTransformed_HistogramReadlength.~{plot_type}"
        File NumberOfReads_Over_Time = "NumberOfReads_Over_Time.~{plot_type}"
        File TimeLengthViolinPlot = "TimeLengthViolinPlot.~{plot_type}"
        File TimeQualityViolinPlot = "TimeQualityViolinPlot.~{plot_type}"
        File TimeSequencingSpeed_ViolinPlot = "TimeSequencingSpeed_ViolinPlot.~{plot_type}"
        File Weighted_HistogramReadlength = "Weighted_HistogramReadlength.~{plot_type}"
        File Weighted_LogTransformed_HistogramReadlength = "Weighted_LogTransformed_HistogramReadlength.~{plot_type}"
        File Yield_By_Length = "Yield_By_Length.~{plot_type}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          "~{num_cpus}",
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "quay.io/biocontainers/nanoplot:1.28.0--py_0"
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
