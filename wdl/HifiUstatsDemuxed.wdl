version 1.0

import "tasks/NanoPlot.wdl"
import "tasks/utils/GeneralUtils.wdl"
import "tasks/Finalize.wdl" as FF

workflow HifiUstatsDemuxed {

    input {
        String smrtcell_data_dir
        String movie

        String gcs_out_root_dir
    }
    String workflow_name = "PreprocessDemultiplexedCCSedSMRTCell"

    call GetAllHifiBams { input: smrtcell_data_dir = smrtcell_data_dir }

    call MergeUBams { input: bams = GetAllHifiBams.bam_paths, prefix = movie }

    call NanoPlot.NanoPlotFromUnAligned as Metrics {input: unaligned_file = MergeUBams.merged_ubam, format='ubam'}

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_name}/~{movie}/"
    call FF.FinalizeToFile as FinalizeMergedHifiBam {
        input:
            file = MergeUBams.merged_ubam,
            outdir = outdir
    }

    call GeneralUtils.TarGZFiles as saveAlnMetrics {
        input:
            files = flatten([[Metrics.stats], Metrics.plots]),
            name = "alignment.metrics"
    }
    call FF.FinalizeToFile as FinalizeNanoplotFiles {
        input:
            file = saveAlnMetrics.you_got_it,
            outdir = outdir
    }

    output {
        File merged_hifi_bam = FinalizeMergedHifiBam.gcs_path
        Map[String, Float] hifi_stats_map = Metrics.stats_map
        File merged_hifi_bam_nanoplot = FinalizeNanoplotFiles.gcs_path
    }
}

task GetAllHifiBams {
    input {
        String smrtcell_data_dir
    }

    command <<<
        set -eux
        gsutil ls -r ~{smrtcell_data_dir} | \
            grep ".bam$" | \
            grep -F 'hifi_reads' \
            > bam.paths.txt
    >>>

    output {
        Array[String] bam_paths = read_lines("bam.paths.txt")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task MergeUBams {
    input {
        Array[File] bams
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bams:   "input array of BAMs to be merged"
        prefix: "[default-valued] prefix for output BAM"
    }

    Int disk_size = 1 + 4*ceil(size(bams, "GB"))

    command <<<
        set -euxo pipefail

        samtools merge \
            -p -c --no-PG \
            -@ 2 \
            -o "~{prefix}.hifi_reads.merged.bam" \
            ~{sep=" " bams}
    >>>

    output {
        File merged_ubam = "~{prefix}.hifi_reads.merged.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
