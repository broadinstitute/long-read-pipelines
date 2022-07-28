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

        ((samtools view -H ~{bam} | grep '^@PG' | grep -w 'PN:ccs') 2>/dev/null) | \
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
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
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

task ShardLongReads {
    input {
        File unaligned_bam
        File unaligned_pbi

        Int num_shards = 300
        Int num_threads = 8

        Boolean drop_per_base_N_pulse_tags = false

        String prefix = "shard"

        String zones = "us-central1-c us-central1-b"
        Int? num_ssds

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        # when running large scale workflows, we sometimes see errors like the following
        #   A resource limit has delayed the operation: generic::resource_exhausted: allocating: selecting resources: selecting region and zone:
        #   no available zones: 2763 LOCAL_SSD_TOTAL_GB (738/30000 available) usage too high
        zones: "select which zone (GCP) to run this task"
    }

    Int disk_size = if defined(num_ssds) then 1 + 375*select_first([num_ssds]) else 1+3*ceil(size(unaligned_bam, "GB") + size(unaligned_pbi, "GB"))
    Int mem = ceil(25*size(unaligned_pbi, "MB")/1000)
    String ex = if drop_per_base_N_pulse_tags then "-x fi,fp,ri,rp" else ""

    command <<<
        set -x

        samtools view -c ~{unaligned_bam} # to check if file is truncated

        python3 /usr/local/bin/shard_bam.py ~{ex} \
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
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.34"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CCS {
    input {
        File subreads

        Boolean all       = true    # see https://ccs.how/faq/mode-all.html for details
        Boolean kinetics  = false   # see https://ccs.how/faq/sqiie.html for details
        Boolean by_strand = false

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(subreads, "GB"))
    String bn = basename(subreads, ".bam")

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        # Move the file from the UUID share to the current folder.
        # This will remove the UUID from the file path and allow call caching to work.
        infile=$( basename ~{subreads} )
        mv ~{subreads} $infile

        # Run CCS:
        ccs ~{true='--all' false='' all} \
            ~{true='--all-kinetics --subread-fallback' false='' kinetics} \
            ~{true='--by-strand' false='' by_strand} \
            --num-threads $num_core \
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.34"
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

task ExtractHifiReads {
    input {
        File bam

        String sample_name
        String library

        String prefix = "hifi"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        sample_name: "we always rely explicitly on input SM name"
        library: "this will override the LB: entry on the @RG line"
    }

    Int disk_size = 1 + 3*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        extracthifi ~{bam} ~{prefix}.tmp.bam

        samtools view --no-PG -H ~{prefix}.tmp.bam > header.txt
        awk '$1 ~ /^@RG/' header.txt > rg_line.txt
        if ! grep -qF "LB:" rg_line.txt; then
            sed -i "s/$/LB:tbd/" rg_line.txt
        fi
        if ! grep -qF "SM:" rg_line.txt; then
            sed -i "s/$/SM:tbd/" rg_line.txt
        fi
        # fix LB:
        awk -v lib="~{library}" -F '\t' 'BEGIN {OFS="\t"} { for (i=1; i<=NF; ++i) { if ($i ~ "LB:") $i="LB:"lib } print}' \
            rg_line.txt \
            > fixed_rg_line.lb.txt
        # fix SM:
        awk -v lib="~{sample_name}" -F '\t' 'BEGIN {OFS="\t"} { for (i=1; i<=NF; ++i) { if ($i ~ "SM:") $i="SM:"lib } print}' \
            fixed_rg_line.lb.txt \
            > fixed_rg_line.txt
        sed -n '/@RG/q;p' header.txt > first_half.txt
        sed -n '/@RG/,$p' header.txt | sed '1d' > second_half.txt

        cat first_half.txt fixed_rg_line.txt second_half.txt > fixed_header.txt

        date
        samtools reheader fixed_header.txt ~{prefix}.tmp.bam > ~{prefix}.bam
        date
    >>>

    output {
        File hifi_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.35"
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

task MergeCCSReports {
    input {
        Array[File] reports
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(reports, "GB"))

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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.35"
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

    Int disk_size = 1 + 2*ceil(size(subreads, "GB") + size(consensus, "GB"))

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

    Int disk_size = 1 + 4*ceil(size(bam, "GB"))

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
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
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

task RefineTranscriptReads {
    input {
        File bam
        File barcode_file
        String prefix = "flnc"
        Boolean require_polya = true

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(bam, "GB"))

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

    Int disk_size = 1 + 4*ceil(size(bam, "GB"))

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

    Int disk_size = 1 + 2*ceil(size([bam, subreads_bam], "GB"))

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
        String? library
        String map_preset

        Boolean drop_per_base_N_pulse_tags

        String prefix = "out"
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        sample_name: "we always rely explicitly on input SM name"
        library: "this will override the LB: entry on the @RG line"
    }

    String median_filter = if map_preset == "SUBREAD" then "--median-filter" else ""
    String extra_options = if drop_per_base_N_pulse_tags then " --strip " else ""

    Boolean fix_library_entry = if defined(library) then true else false

    Int disk_size = 1 + 10*ceil(size(bam, "GB") + size(ref_fasta, "GB"))
    Int cpus = 16
    Int mem = 96

    command <<<
        set -euxo pipefail

        pbmm2 align ~{bam} ~{ref_fasta} ~{prefix}.pre.bam \
            --preset ~{map_preset} \
            ~{median_filter} \
            --sample ~{sample_name} \
            ~{extra_options} \
            --sort \
            --unmapped

        if ~{fix_library_entry}; then
            mv ~{prefix}.pre.bam ~{prefix}.pre.tmp.bam
            samtools view --no-PG -H ~{prefix}.pre.tmp.bam > header.txt
            awk '$1 ~ /^@RG/' header.txt > rg_line.txt
            awk -v lib="~{library}" -F '\t' 'BEGIN {OFS="\t"} { for (i=1; i<=NF; ++i) { if ($i ~ "LB:") $i="LB:"lib } print}' \
                rg_line.txt \
                > fixed_rg_line.txt
            sed -n '/@RG/q;p' header.txt > first_half.txt
            sed -n '/@RG/,$p' header.txt | sed '1d' > second_half.txt

            cat first_half.txt fixed_rg_line.txt second_half.txt > fixed_header.txt

            date
            samtools reheader fixed_header.txt ~{prefix}.pre.tmp.bam > ~{prefix}.pre.bam
            rm ~{prefix}.pre.tmp.bam
            date
        fi

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

    Int disk_size = 1 + 4*ceil(size(bam, "GB"))

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

    Int disk_size = 1 + 2*ceil(size(report, "GB"))

    command <<<
        set -euxo pipefail

        cat ~{report} | grep 'ZMWs input' | awk -F": " '{ print $2 }' > zmws_input.txt
        cat ~{report} | grep 'ZMWs pass filters' | awk -F": " '{ print $2 }' | awk '{ print $1 }' > zmws_pass_filters.txt
        cat ~{report} | grep 'ZMWs fail filters' | awk -F": " '{ print $2 }' | awk '{ print $1 }' > zmws_fail_filters.txt
        cat ~{report} | grep 'ZMWs shortcut filters' | awk -F": " '{ print $2 }' | awk '{ print $1 }' > zmws_shortcut_filters.txt
        cat ~{report} | grep 'ZMWs pass filters' | awk -F": " '{ print $2 }' | awk '{ print $2 }' | sed 's/[()%]//g' > zmws_pass_filters_pct.txt
        cat ~{report} | grep 'ZMWs fail filters' | awk -F": " '{ print $2 }' | awk '{ print $2 }' | sed 's/[()%]//g' > zmws_fail_filters_pct.txt
        cat ~{report} | grep 'ZMWs shortcut filters' | awk -F": " '{ print $2 }' | awk '{ print $2 }' | sed 's/[()%]//g' > zmws_shortcut_filters_pct.txt
    >>>

    output {
        Float zmws_input = read_float("zmws_input.txt")
        Float zmws_pass_filters = read_float("zmws_pass_filters.txt")
        Float zmws_fail_filters = read_float("zmws_fail_filters.txt")
        Float zmws_shortcut_filters = read_float("zmws_shortcut_filters.txt")
        Float zmws_pass_filters_pct = read_float("zmws_pass_filters_pct.txt")
        Float zmws_fail_filters_pct = read_float("zmws_fail_filters_pct.txt")
        Float zmws_shortcut_filters_pct = read_float("zmws_shortcut_filters_pct.txt")
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

    Int disk_size = 1 + 2*ceil(size(xml, "GB"))

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

    Int disk_size = 1 + 2*ceil(size(pbi, "GB"))

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
        mem_gb:             16,
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
