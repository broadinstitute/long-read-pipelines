version 1.0

import "Structs.wdl"

task FindBams {
    input {
        String gcs_input_dir

        RuntimeAttr? runtime_attr_override
    }

    String indir = sub(gcs_input_dir, "/$", "")

    command <<<
        # Commenting out this line to make it able to run even if there are no subreads / reads files:
        # set -euxo pipefail

        # Get our file lists here:
        gsutil ls ~{indir}/**subreads.bam | sort > subread_bams.txt
        gsutil ls ~{indir}/**.reads.bam | sort > ccs_bams.txt

        # Check our file lists to see which are populated:
        if [[ $( wc -l subread_bams.txt | awk '{print $1}' ) -gt 0 ]] ; then
            echo 'Found subreads!'
            echo "true" > has_subreads.txt
        else
            echo 'Did not find subreads!'
            echo "false" > has_subreads.txt
        fi

        if [[ $( wc -l ccs_bams.txt | awk '{print $1}' ) -gt 0 ]] ; then
            echo 'Found CCS reads!'
            echo "true" > has_ccs_reads.txt
        else
            echo 'Did not find CCS reads!'
            echo "false" > has_ccs_reads.txt
        fi
    >>>

    output {
        Array[String] subread_bams = read_lines("subread_bams.txt")
        Array[String] ccs_bams = read_lines("ccs_bams.txt")

        Boolean has_subreads = read_boolean("has_subreads.txt")
        Boolean has_ccs_reads = read_boolean("has_ccs_reads.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
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

task FindZmwStatsJsonGz {
    input {
        String gcs_input_dir

        RuntimeAttr? runtime_attr_override
    }

    String indir = sub(gcs_input_dir, "/$", "")

    command <<<
        # Commenting out this line to make it able to run even if there are no subreads / reads files:
        # set -euxo pipefail

        # Get our file lists here:
        gsutil ls ~{indir}/**.zmw_metrics.json.gz | sort > files.txt
    >>>

    output {
        File zmw_stats_json_gz = read_lines("files.txt")[0]
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
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

task FindCCSReport {
    input {
        String gcs_input_dir

        RuntimeAttr? runtime_attr_override
    }

    String indir = sub(gcs_input_dir, "/$", "")

    command <<<
        set -euxo pipefail
        gsutil ls ~{indir}/**.ccs_reports.txt | sort > ccs_reports.txt
    >>>

    output {
        # There should be only one CCS report here:
        File ccs_report = read_lines("ccs_reports.txt")[0]
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
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
        String subread_bam

        String? bam_suffix

        RuntimeAttr? runtime_attr_override
    }

    String gcs_dir = sub(subread_bam, basename(subread_bam), "")

    String bam_suffix_arg = if defined(bam_suffix) then " --BS " else ""

    command <<<
        set -x

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

        # We need to update detect_run_info.py to make it sanitize fields.
        # The `sed` statement here is a hack to get around an issue.
        python /usr/local/bin/detect_run_info.py ~{gcs_dir} ~{bam_suffix_arg}~{default="" sep=" --BS " bam_suffix} > run_info.txt
    >>>

    output {
        Map[String, String] run_info = read_map("run_info.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.30"
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

task RemoveKineticsTags {
    input {
        File bam
        String prefix = "kinetics_removed"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(bam, "GB"))
    command <<<
        set -x

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        samtools view -hb -@$np -x fi -x ri -x fp -x rp -x fn -x rn ~{bam} > ~{prefix}.bam
    >>>

    output {
        File bam_file = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.21"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.21"
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

        Int min_passes = 3
        Float min_snr = 2.5
        Int min_length = 10
        Int max_length = 50000
        Float min_rq = 0.99
        Boolean by_strand = false

        Int disk_space_scale_factor = 2

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = disk_space_scale_factor*ceil(size(subreads, "GB"))

    command <<<
        set -euxo pipefail

        # Move the file from the UUID share to the current folder.
        # This will remove the UUID from the file path and allow call caching to work.
        infile=$( basename ~{subreads} )
        mv ~{subreads} $infile

        # Get the number of processors we're runnign with:
        num_procs=$(cat /proc/cpuinfo | grep processor | awk '{print $NF}' | sort -n | tail -n1)
        let num_procs=$num_procs+1

        # Run CCS:
        ccs --min-passes ~{min_passes} \
            --min-snr ~{min_snr} \
            --min-length ~{min_length} \
            --max-length ~{max_length} \
            --min-rq ~{min_rq} \
            --num-threads $num_procs \
            --report-file ccs_report.txt \
            ~{if by_strand then "--by-strand" else ""} $infile ccs_unmapped.bam
    >>>

    output {
        File consensus = "ccs_unmapped.bam"
        File report = "ccs_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.21"
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

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(reports, "GB")) + 1

    command <<<
        set -euxo pipefail

        python /usr/local/bin/merge_ccs_reports.py ~{sep=' ' reports} > ccs_report.txt
    >>>

    output {
        File report = "ccs_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.21"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.22"
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

task CountZMWs {
    input {
        File bam
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bam, "GB")) + 1

    command <<<
        samtools view ~{bam} | awk 'BEGIN{FS="/";OFS="\t"}{print $2}'  | uniq | sort -n | uniq | wc -l
    >>>

    output {
        Int num_zmws = read_int(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.21"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.21"
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

    Int disk_size = 1 + 2*ceil(size(report, "GB"))

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
        mem_gb:             8,
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.21"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.21"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.21"
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

    Int disk_size = 1 + 10*ceil(size(bam, "GB") + size(ref_fasta, "GB"))

    command <<<
        set -euxo pipefail

        pbmm2 align ~{bam} ~{ref_fasta} ~{prefix}.bam --preset ~{map_preset} --sort

        find . -type f -exec ls -lah {} \;
    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.21"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.21"
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

task CollectZmwSubreadStats {
    input {
        File subreads
        String prefix = "zmw_subread_stats"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(subreads, "GB"))

    command <<<
        set -euxo pipefail

        python3 /usr/local/bin/collect_zmw_subread_stats.py -o ~{prefix}.tsv -b ~{subreads}
    >>>

    output {
        File zmw_subread_stats = "~{prefix}.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.24.1"
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

task CollectPolymeraseReadLengths {
    input {
        String gcs_input_dir
        File subreads
        String prefix = "polymerase_read_lengths"

        RuntimeAttr? runtime_attr_override
    }

    String indir = sub(gcs_input_dir, "/$", "")
    Int disk_size = 2*ceil(size(subreads, "GB"))

    command <<<
        set -euxo pipefail

        # Find our corresponding scraps file:
        subreads_base_name=$( basename ~{subreads} | sed 's#.subreads.bam$##' )
        gsutil ls ~{gcs_input_dir}/**scraps.bam | grep "$subreads_base_name" > scraps_file_name.txt
        if [[ $( wc -l scraps_file_name.txt | awk '{print $1}' ) -ne 1 ]] ; then
            echo "ERROR: Cannot match scraps file to given subreads file: $subreads_base_name" 1>&2
            echo "       Multiple files returned:" 1>&2
            cat scraps_file_name.txt 1>&2
            return 1
        fi

        # Localize the scraps so we can read them:
        remote_scraps="$(head -n1 scraps_file_name.txt)"
        local_scraps="$(basename $remote_scraps)"
        echo "Localizing scraps: $remote_scraps -> $local_scraps"
        gsutil -m cp $remote_scraps .

        # Calculate our polymerase read lengths:
        python3 /usr/local/bin/collect_polymerase_read_lengths.py -o ~{prefix}.tsv -s ~{subreads} -p $local_scraps
    >>>

    output {
        File polymerase_read_lengths_tsv = "~{prefix}.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.24.1"
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

    Int disk_size = 2*ceil(size(bam, "GB"))

    String base_name = basename(bam)

    command <<<
        set -euxo pipefail

        # Run PBIndex:
        pbindex ~{bam}
        mv ~{bam}.pbi ~{base_name}.pbi
    >>>

    output {
        File pbindex = "~{base_name}.pbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,             # 1 preemptible try here so that we have even money with our local SSDs below.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.21"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"  # LOCAL here is a local SSD - much faster and even money with normal disk if preemptible
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CopyMetadataFilesToNewDir {
    input {
        String input_gs_path
        String dest_gs_path

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        gsutil -m cp ~{input_gs_path}/*.xml ~{input_gs_path}/*.log ~{input_gs_path}/*.json ~{input_gs_path}/*.txt  ~{input_gs_path}/*.transferdone ~{input_gs_path}/*.json.gz ~{dest_gs_path}/
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
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