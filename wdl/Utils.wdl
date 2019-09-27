version 1.0

import "Structs.wdl"

# TODO: describe purpose
task DetectRunInfo {
    input {
        String gcs_dir
        String? sample_name

        RuntimeAttr? runtime_attr_override
    }

    String SM = if defined(sample_name) then "--SM " + sample_name else ""

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        python /usr/local/bin/detect_run_info.py ~{SM} ~{gcs_dir} > run_info.txt
        gsutil ls ~{gcs_dir} | grep -v scraps | grep -e '\.bam$' -e '\.f\(ast\)\?q\(\.gz\)\?$' > files.txt
    >>>

    output {
        File run_info_file = "run_info.txt"
        File fofn = "files.txt"
        Map[String, String] run_info = read_map("run_info.txt")
        Array[String] files = read_lines("files.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1, 
        mem_gb:             1, 
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "kgarimella/lr-align:0.01.17"
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

# TODO: describe purpose
task PrepareRun {
    input {
        Array[File] files
        
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(files, "GB"))

    command <<<
        set -euxo pipefail

        python /usr/local/bin/prepare_run.py ~{sep=' ' files}
    >>>

    output {
        File unmapped_bam = "unmapped.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2, 
        mem_gb:             2, 
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        3,
        docker:             "kgarimella/lr-align:0.01.17"
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

# TODO: investigate if samtools is better for this task
# Sort BAM file by coordinate order
# copied from "dsde_pipelines_tasks/BamProcessing.wdl", with
#    customization on the runtime block, and "preemptible_tries" taken out
task SortSam {
    input {
      File input_bam
      String output_bam_basename
      Int compression_level

      RuntimeAttr? runtime_attr_override
    }
    # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
    # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
    Float sort_sam_disk_multiplier = 3.25
    Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20

    command {
      java -Dsamjdk.compression_level=~{compression_level} -Xms4000m -jar /usr/gitc/picard.jar \
        SortSam \
        INPUT=~{input_bam} \
        OUTPUT=~{output_bam_basename}.bam \
        SORT_ORDER="coordinate" \
        CREATE_INDEX=true \
        CREATE_MD5_FILE=true \
        MAX_RECORDS_IN_RAM=300000
    }
    output {
      File output_bam = "~{output_bam_basename}.bam"
      File output_bam_index = "~{output_bam_basename}.bai"
      File output_bam_md5 = "~{output_bam_basename}.bam.md5"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             5,
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
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