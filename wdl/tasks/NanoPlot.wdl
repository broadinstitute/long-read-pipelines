version 1.0

import "Structs.wdl"

task NanoPlotFromSummary {
    input {
        Array[File] summary_files

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(summary_files, "GB"))

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        NanoPlot -t ${num_core} \
                 -c orangered \
                 --N50 \
                 --summary "~{sep=' ' summary_files}"
    >>>

    output {
        File stats = "NanoStats.txt"

        Array[File] plots = glob("*.png")
        File ActivePores_Over_Time = "ActivePores_Over_Time.png"
        File ActivityMap_ReadsPerChannel = "ActivityMap_ReadsPerChannel.png"
        File CumulativeYieldPlot_Gigabases = "CumulativeYieldPlot_Gigabases.png"
        File CumulativeYieldPlot_NumberOfReads = "CumulativeYieldPlot_NumberOfReads.png"
        File LengthvsQualityScatterPlot_dot = "LengthvsQualityScatterPlot_dot.png"
        File LengthvsQualityScatterPlot_kde = "LengthvsQualityScatterPlot_kde.png"
        File Non_weightedHistogramReadlength = "Non_weightedHistogramReadlength.png"
        File Non_weightedLogTransformed_HistogramReadlength = "Non_weightedLogTransformed_HistogramReadlength.png"
        File NumberOfReads_Over_Time = "NumberOfReads_Over_Time.png"
        File TimeLengthViolinPlot = "TimeLengthViolinPlot.png"
        File TimeQualityViolinPlot = "TimeQualityViolinPlot.png"
        File TimeSequencingSpeed_ViolinPlot = "TimeSequencingSpeed_ViolinPlot.png"
        File WeightedHistogramReadlength = "WeightedHistogramReadlength.png"
        File WeightedLogTransformed_HistogramReadlength = "WeightedLogTransformed_HistogramReadlength.png"
        File Yield_By_Length = "Yield_By_Length.png"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/biocontainers/nanoplot:1.35.5--pyhdfd78af_0"
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
