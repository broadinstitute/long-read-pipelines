version 1.0

import "../../structs/Structs.wdl"
import "../../tasks/Utility/Finalize.wdl" as FF

workflow UnalignedMetrics {
    input {
        File unaligned_bam

        String per
        String type
        String label

        String? gcs_output_dir
    }

    call FlagStats as UnalignedFlagStats {
        input:
            bam = unaligned_bam,
    }

    call ReadMetrics as UnalignedReadMetrics {
        input:
            bam = unaligned_bam,
    }

    call ReadNamesAndLengths {
        input:
            bam = unaligned_bam
    }

    if (defined(gcs_output_dir)) {
        String outdir = sub(sub(gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://") + "/metrics/~{per}_~{type}_~{label}"

        call FF.FinalizeToDir as FFYieldUnaligned {
            input:
                outdir = outdir + "/yield_unaligned/",
                files = [
                    UnalignedFlagStats.flag_stats,
                    UnalignedReadMetrics.np_hist,
                    UnalignedReadMetrics.range_gap_hist,
                    UnalignedReadMetrics.zmw_hist,
                    UnalignedReadMetrics.prl_counts,
                    UnalignedReadMetrics.prl_hist,
                    UnalignedReadMetrics.prl_nx,
                    UnalignedReadMetrics.prl_yield_hist,
                    UnalignedReadMetrics.rl_counts,
                    UnalignedReadMetrics.rl_hist,
                    UnalignedReadMetrics.rl_nx,
                    UnalignedReadMetrics.rl_yield_hist
                ]
        }

        call FF.FinalizeToDir as FFReadNamesAndLengths { input: outdir = outdir + "/read_names_and_lengths/", files = [ ReadNamesAndLengths.read_names_and_lengths ] }
    }

    output {
        File unaligned_flag_stats = UnalignedFlagStats.flag_stats

        File unaligned_np_hist = UnalignedReadMetrics.np_hist
        File unaligned_range_gap_hist = UnalignedReadMetrics.range_gap_hist
        File unaligned_zmw_hist = UnalignedReadMetrics.zmw_hist
        File unaligned_prl_counts = UnalignedReadMetrics.prl_counts
        File unaligned_prl_hist = UnalignedReadMetrics.prl_hist
        File unaligned_prl_nx = UnalignedReadMetrics.prl_nx
        File unaligned_prl_yield_hist = UnalignedReadMetrics.prl_yield_hist
        File unaligned_rl_counts = UnalignedReadMetrics.rl_counts
        File unaligned_rl_hist = UnalignedReadMetrics.rl_hist
        File unaligned_rl_nx = UnalignedReadMetrics.rl_nx
        File unaligned_rl_yield_hist = UnalignedReadMetrics.rl_yield_hist
    }
}

task FlagStats {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools flagstat ~{bam} > ~{basename}.flag_stats.txt
    >>>

    output {
        File flag_stats = "~{basename}.flag_stats.txt"
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

task ReadMetrics {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        java -jar /usr/local/bin/gatk.jar ComputeLongReadMetrics -I ~{bam} -O ~{basename}.read_metrics -DF WellformedReadFilter
    >>>

    output {
        File np_hist = "~{basename}.read_metrics.np_hist.txt"
        File range_gap_hist = "~{basename}.read_metrics.range_gap_hist.txt"
        File zmw_hist = "~{basename}.read_metrics.zmw_hist.txt"
        File prl_counts = "~{basename}.read_metrics.prl_counts.txt"
        File prl_hist = "~{basename}.read_metrics.prl_hist.txt"
        File prl_nx = "~{basename}.read_metrics.prl_nx.txt"
        File prl_yield_hist = "~{basename}.read_metrics.prl_yield_hist.txt"
        File rl_counts = "~{basename}.read_metrics.rl_counts.txt"
        File rl_hist = "~{basename}.read_metrics.rl_hist.txt"
        File rl_nx = "~{basename}.read_metrics.rl_nx.txt"
        File rl_yield_hist = "~{basename}.read_metrics.rl_yield_hist.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             40,
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

task ReadNamesAndLengths {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools view ~{bam} | awk '{ print $1, length($10) }' | gzip -1 > ~{basename}.read_names_and_lengths.txt.gz
    >>>

    output {
        File read_names_and_lengths = "~{basename}.read_names_and_lengths.txt.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
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

