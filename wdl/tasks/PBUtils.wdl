version 1.0

import "Structs.wdl"

task FindBams {
    input {
        String gcs_input_dir

        RuntimeAttr? runtime_attr_override
    }

    String indir = sub(gcs_input_dir, "/$", "")

    command <<<
        set -euxo pipefail

        gsutil ls "~{indir}/**subreads.bam" > subread_bams.txt
    >>>

    output {
        Array[String] subread_bams = read_lines("subread_bams.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
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

task GetRunInfo {
    input {
        String bam
        String SM

        RuntimeAttr? runtime_attr_override
    }

    String gcs_dir = sub(bam, basename(bam), "")

    command <<<
        set -x

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        python /usr/local/bin/detect_run_info.py --SM ~{SM} ~{gcs_dir} > run_info.txt

        ((samtools view -H ~{bam} | grep '^@PG[[:space:]]\+ID:ccs-') 2>/dev/null) | \
            wc -l | \
            sed 's/0/false/' | \
            sed 's/1/true/' \
            > status.txt
    >>>

    output {
        Map[String, String] run_info = read_map("run_info.txt")
        Boolean is_corrected = read_boolean("status.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.29"
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

task ShardLongReads {
    input {
        File unaligned_bam
        File unaligned_pbi

        Int num_shards = 300
        Int num_threads = 8

        String prefix = "shard"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3*ceil(size(unaligned_bam, "GB") + size(unaligned_pbi, "GB"))
    Int mem = ceil(25*size(unaligned_pbi, "MB")/1000)

    command <<<
        set -x

        python3 /usr/local/bin/shard_bam.py \
            -n ~{num_shards} \
            -t ~{num_threads} \
            -i ~{unaligned_pbi} \
            -p ~{prefix} \
            ~{unaligned_bam}
    >>>

    output {
        Array[File] unmapped_shards = glob("~{prefix}*.bam")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_threads,
        mem_gb:             mem,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.29"
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

task CCS {
    input {
        File subreads

        Boolean by_strand = false
        Int cpus = 4

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(subreads, "GB"))
    String bn = basename(subreads, ".bam")

    command <<<
        set -euxo pipefail

        # Move the file from the UUID share to the current folder.
        # This will remove the UUID from the file path and allow call caching to work.
        infile=$( basename ~{subreads} )
        mv ~{subreads} $infile

        # Run CCS:
        ccs --all --all-kinetics --subread-fallback ~{if by_strand then "--by-strand" else ""} \
            --num-threads ~{cpus} \
            --log-level INFO \
            --log-file ~{bn}.ccs.log \
            --stderr-json-log \
            --suppress-reports \
            --report-file ~{bn}.ccs_reports.txt \
            --report-json ~{bn}.ccs_reports.json \
            --metrics-json ~{bn}.zmw_metrics.json.gz \
            --hifi-summary-json ~{bn}.hifi_summary.json \
            $infile ~{bn}.ccs_unmapped.bam
    >>>

    output {
        File consensus = "~{bn}.ccs_unmapped.bam"
        File report = "~{bn}.ccs_reports.txt"
        File report_json = "~{bn}.ccs_reports.json"
        File metrics_json = "~{bn}.zmw_metrics.json.gz"
        File hifi_summary_json = "~{bn}.hifi_summary.json"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             12,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.29"
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

task MergeCCSReports {
    input {
        Array[File] reports
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(reports, "GB")) + 1

    command <<<
        set -euxo pipefail

        python /usr/local/bin/merge_ccs_reports.py ~{sep=' ' reports} > ~{prefix}.ccs_report.txt
    >>>

    output {
        File report = "~{prefix}.ccs_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.29"
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

task ExtractUncorrectedReads {
    input {
        File subreads
        File consensus

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(subreads, "GB") + size(consensus, "GB"))

    command <<<
        set -euxo pipefail

        python3 /usr/local/bin/extract_uncorrected_reads.py -o ~{prefix}.bam ~{subreads} ~{consensus}
    >>>

    output {
        File uncorrected = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.29"
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

task Demultiplex {
    input {
        File bam
        File barcode_file
        String prefix           = "demux"
        Boolean ccs             = false
        Boolean isoseq          = false
        Boolean peek_guess      = false
        Boolean dump_removed    = false
        Boolean split_bam_named = false
        Int peek                = 0
        Int min_score           = 0
        Int guess               = 0
        Int guess_min_count     = 0

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bam, "GB")) + 1

    command <<<
        set -euxo pipefail

        lima \
            ~{if ccs then "--ccs" else ""} \
            ~{if isoseq then "--isoseq" else ""} \
            ~{if peek_guess then "--peek-guess" else ""} \
            ~{if guess > 0 then "--guess ~{guess}" else ""} \
            ~{if guess_min_count > 0 then "--guess-min-count ~{guess_min_count}" else ""} \
            ~{if peek > 0 then "--peek ~{peek}" else ""} \
            ~{if dump_removed then "--dump-removed" else ""} \
            ~{if split_bam_named then "--split-bam-named" else ""} \
            ~{bam} \
            ~{barcode_file} \
            ~{prefix}.bam

        find . -type f -exec ls -lah {} \;
    >>>

    output {
        Array[File] demux_bams = glob("~{prefix}.*.bam")
        File counts = "~{prefix}.lima.counts"
        #File guess = "~{prefix}.lima.guess"
        File report = "~{prefix}.lima.report"
        File summary = "~{prefix}.lima.summary"
        File? clips = "~{prefix}.lima.clips"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.29"
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

task MakeDetailedDemultiplexingReport {
    input {
        File report
        String type = "png"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(report, "GB"))

    command <<<
        set -euxo pipefail

        Rscript /lima_report_detail.R ~{report} ~{type}
    >>>

    output {
        Array[File] report_files = glob("detail_*~{type}")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
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

task MakeSummarizedDemultiplexingReport {
    input {
        File report

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(report, "GB"))

    command <<<
        set -euxo pipefail

        Rscript /lima_report_summary.R ~{report}
    >>>

    output {
        Array[File] report_files = glob("summary_*.png")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
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

task MakePerBarcodeDemultiplexingReports {
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
        boot_disk_gb:       10,
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

task RefineTranscriptReads {
    input {
        File bam
        File barcode_file
        String prefix = "flnc"
        Boolean require_polya = true

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bam, "GB")) + 1

    command <<<
        set -euxo pipefail

        isoseq3 refine ~{bam} ~{barcode_file} ~{prefix}.bam ~{true='--require-polya' false='' require_polya}
    >>>

    output {
        File refined_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.29"
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

task ClusterTranscripts {
    input {
        File bam
        String prefix = "clustered"
        Boolean use_qvs = true

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bam, "GB")) + 1

    command <<<
        set -euxo pipefail

        isoseq3 cluster ~{bam} ~{prefix}.bam --verbose ~{true='--use-qvs' false='' use_qvs}
    >>>

    output {
        File clustered_bam = "~{prefix}.bam"
        File clustered_pbi = "~{prefix}.bam.pbi"
        File hq_fasta = "~{prefix}.hq.fasta.gz"
        File hq_bam = "~{prefix}.hq.bam"
        File hq_pbi = "~{prefix}.hq.bam.pbi"
        File lq_fasta = "~{prefix}.lq.fasta.gz"
        File lq_bam = "~{prefix}.lq.bam"
        File lq_pbi = "~{prefix}.lq.bam.pbi"
        File cluster = "~{prefix}.cluster"
        File cluster_report_csv = "~{prefix}.cluster_report.csv"
        File transcriptset_xml = "~{prefix}.transcriptset.xml"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          64,
        mem_gb:             70,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.29"
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

task PolishTranscripts {
    input {
        File bam
        File subreads_bam
        File subreads_pbi
        String prefix = "polished"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([bam, subreads_bam], "GB"))

    command <<<
        set -euxo pipefail

        isoseq3 polish ~{bam} ~{subreads_bam} ~{prefix}.bam
    >>>

    output {
        File polished_bam = "~{prefix}.bam"
        File polished_fastq = "~{prefix}.hq.fastq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          24,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.29"
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

task Align {
    input {
        File bam
        File ref_fasta

        String sample_name
        String map_preset

        String prefix = "out"
        RuntimeAttr? runtime_attr_override
    }

    String median_filter = if map_preset == "SUBREAD" then "--median-filter" else ""

    Int disk_size = 1 + 10*ceil(size(bam, "GB") + size(ref_fasta, "GB"))
    Int cpus = 4
    Int mem = 30

    command <<<
        set -euxo pipefail

        pbmm2 align ~{bam} ~{ref_fasta} ~{prefix}.pre.bam --preset ~{map_preset} ~{median_filter} --sort

        samtools calmd -b --no-PG ~{prefix}.pre.bam ~{ref_fasta} > ~{prefix}.bam
        samtools index ~{prefix}.bam
    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             mem,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.29"
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

task PBIndex {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        mv ~{bam} ~{basename(bam)}

        pbindex ~{basename(bam)}
    >>>

    output {
        File pbi = "~{basename(bam)}.pbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.29"
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

task CollapseTranscripts {
    input {
        File bam
        String prefix = "out"
        Boolean use_qvs = true

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bam, "GB")) + 1

    command <<<
        set -euxo pipefail

        isoseq3 collapse ~{bam} ~{prefix}.gff
    >>>

    output {
        File gff = "~{prefix}.gff"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.29"
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

task SummarizeCCSReport {
    input {
        File report

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(report, "GB"))

    command <<<
        set -euxo pipefail

        cat ~{report} | grep 'ZMWs input' | awk -F": " '{ print $2 }' > zmws_input.txt
        cat ~{report} | grep 'ZMWs pass filters' | awk -F": " '{ print $2 }' | awk '{ print $1 }' > zmws_pass_filters.txt
        cat ~{report} | grep 'ZMWs fail filters' | awk -F": " '{ print $2 }' | awk '{ print $1 }' > zmws_fail_filters.txt
        cat ~{report} | grep 'ZMWs pass filters' | awk -F": " '{ print $2 }' | awk '{ print $2 }' | sed 's/[()%]//g' > zmws_pass_filters_pct.txt
        cat ~{report} | grep 'ZMWs fail filters' | awk -F": " '{ print $2 }' | awk '{ print $2 }' | sed 's/[()%]//g' > zmws_fail_filters_pct.txt
    >>>

    output {
        Float zmws_input = read_float("zmws_input.txt")
        Float zmws_pass_filters = read_float("zmws_pass_filters.txt")
        Float zmws_fail_filters = read_float("zmws_fail_filters.txt")
        Float zmws_pass_filters_pct = read_float("zmws_pass_filters_pct.txt")
        Float zmws_fail_filters_pct = read_float("zmws_fail_filters_pct.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
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

task SummarizeXMLMetadata {
    input {
        File xml

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(xml, "GB"))

    command <<<
        set -euxo pipefail

        cat ~{xml} | grep '<pbds:TotalLength>' | sed 's/<pbds:TotalLength>//g' | sed 's/<\/pbds:TotalLength>//' | sed 's/\s*//g' > xml_total_length.txt
        cat ~{xml} | grep '<pbds:NumRecords>' | sed 's/<pbds:NumRecords>//g' | sed 's/<\/pbds:NumRecords>//' | sed 's/\s*//g' > xml_num_records.txt
    >>>

    output {
        Float xml_total_length = read_float("xml_total_length.txt")
        Float xml_num_records = read_float("xml_num_records.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
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

task SummarizePBI {
    input {
        File pbi
        Int qual_threshold = 0

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(pbi, "GB"))

    command <<<
        set -euxo pipefail

        python3 /usr/local/bin/compute_pbi_stats.py -q ~{qual_threshold} ~{pbi} | tee map.txt
    >>>

    output {
        Map[String, Float] results = read_map("map.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             12,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.29"
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
