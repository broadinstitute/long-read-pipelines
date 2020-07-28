version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.23/wdl/tasks/Structs.wdl"

task FindBams {
    input {
        String gcs_input_dir

        RuntimeAttr? runtime_attr_override
    }

    String indir = sub(gcs_input_dir, "/$", "")

    command <<<
        set -euxo pipefail

        gsutil ls ~{indir}/**subreads.bam > subread_bams.txt
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

        RuntimeAttr? runtime_attr_override
    }

    String gcs_dir = sub(subread_bam, basename(subread_bam), "")

    command <<<
        set -x

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        python /usr/local/bin/detect_run_info.py ~{gcs_dir} > run_info.txt
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.5"
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

        Int cpus = 4

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(subreads, "GB"))

    command <<<
        set -euxo pipefail

        # Move the file from the UUID share to the current folder.
        # This will remove the UUID from the file path and allow call caching to work.
        infile=$( basename ~{subreads} )
        mv ~{subreads} $infile

        # Run CCS:
        ccs --min-passes ~{min_passes} \
            --min-snr ~{min_snr} \
            --min-length ~{min_length} \
            --max-length ~{max_length} \
            --min-rq ~{min_rq} \
            --num-threads ~{cpus} \
            $infile ccs_unmapped.bam
    >>>

    output {
        File consensus = "ccs_unmapped.bam"
        File report = "ccs_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             40,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.5"
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

task CCSWithClasses {
    input {
        File subreads

        Int min_passes = 3
        Float min_snr = 2.5
        Int min_length = 10
        Int max_length = 50000
        Float min_rq = 0.99

        Int cpus = 4

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(subreads, "GB"))

    command <<<
        set -x

        mkdir shards corrected reports

        java -Dsamjdk.compression_level=0 -jar /usr/local/bin/gatk.jar ShardLongReads -I ~{subreads} -nr 1 -O shards/ -DF WellformedReadFilter --use-jdk-deflater --use-jdk-inflater

        for b in shards/*.bam
        do
            name=`samtools view $b | head -1 | awk -F"/" '{ print $1 "__" $2 }'`

            ccs --min-passes ~{min_passes} \
                --min-snr ~{min_snr} \
                --min-length ~{min_length} \
                --max-length ~{max_length} \
                --min-rq ~{min_rq} \
                --num-threads ~{cpus} \
                --report-file reports/$name.ccs_report.txt \
                $b corrected/ccs.$name.bam
        done

        grep -H '100' reports/*report.txt | grep -v -e 'ZMWs input' -e 'ZMWs filtered' | sed 's/.ccs_report.txt:/\t/' | sed 's/:.*//' | sed 's/ (B)//' | sed 's/reports\///' | sed 's/__/\//' | gzip -1 > ccs_classes.txt.gz
        samtools merge ccs_unmapped.bam corrected/*.bam
        python /usr/local/bin/merge_ccs_reports.py reports/*.txt > ccs_report.txt
    >>>

    output {
        File consensus = "ccs_unmapped.bam"
        File report = "ccs_report.txt"
        File classes = "ccs_classes.txt.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             40,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.5"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.5"
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

task MergeCCSClasses {
    input {
        Array[File] classes

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(classes, "GB")) + 1

    command <<<
        set -euxo pipefail

        cat ~{sep=' ' classes} > ccs_classes.txt.gz
    >>>

    output {
        File classes = "ccs_classes.txt.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.5"
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
