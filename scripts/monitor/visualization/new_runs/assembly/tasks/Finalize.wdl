version 1.0

import "Structs.wdl"

task FinalizeToFile {
    input {
        File file
        String outdir
        String? name
        Boolean gzip_compress = false

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        file: {
            description: "file to finalize",
            localization_optional: true
        }
        outdir: "directory to which files should be uploaded"
        name:   "name to set for uploaded file"
    }

    String base = basename(file)
    String gcs_output_dir = sub(outdir, "/+$", "")
    # THIS IS ABSOLUTELY CRITICAL: DON'T CHANGE TYPE TO FILE, AS CROMWELL WILL TRY TO LOCALIZE THIS NON-EXISTENT FILE
    String gcs_output_file = gcs_output_dir + "/" + select_first([name, base]) + if gzip_compress then ".gz" else ""

    Int disk_size = 2 * ceil(size(file, "GB"))

    command <<<
        set -euxo pipefail

        if ~{gzip_compress}; then
            gsutil cp "~{file}" . && gzip -vk "~{base}" ## need to locatize explicitly
            gsutil -m cp "~{base}.gz" "~{gcs_output_file}"
        else
            gsutil -m cp "~{file}" "~{gcs_output_file}"
        fi

    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.3"
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

task CompressAndFinalize {
    input {
        File file
        String outdir
        String? name
        Boolean gzip_compress = false

        RuntimeAttr? runtime_attr_override
    }

    String base = basename(file)
    String gcs_output_dir = sub(outdir, "/+$", "")
    # THIS IS ABSOLUTELY CRITICAL: DON'T CHANGE TYPE TO FILE, AS CROMWELL WILL TRY TO LOCALIZE THIS NON-EXISTENT FILE
    String gcs_output_file = gcs_output_dir + "/" + select_first([name, base]) + if gzip_compress then ".gz" else ""

    Int disk_size = 2 * ceil(size(file, "GB"))

    command <<<
        set -euxo pipefail
        
        ## need to localize explicitly
        gzip -vkc ~{file} > "~{base}.gz" 
        find ./ -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g'
        gsutil cp "~{base}.gz" "~{gcs_output_file}"

    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.3"
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

task FinalizeAndCompress {
    input {
        Array[File] files
        String outdir

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    String gcs_output_file = sub(outdir, "/+$", "") + "/" + prefix

    Int disk_size = 5 * ceil(size(files, "GB"))

    command <<<
        set -euxo pipefail

        for ff in ~{sep=' ' files};
        do
            base=$(echo "${ff}" | awk -F '/' '{print $NF}')
            gsutil cp "${ff}" . && gzip -vk "${base}" ## need to locatize explicitly
        done
        find ./ -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g' ## emulate tree command
        ls /cromwell_root/*.gz | gsutil -m cp -I "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             7,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.3"
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

task FinalizeToDir {
    input {
        Array[File] files
        String outdir

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        files: {
            description: "files to finalize",
            localization_optional: true
        }
        outdir: "directory to which files should be uploaded"
    }

    String gcs_output_dir = sub(outdir, "/+$", "")

    command <<<
        set -euxo pipefail

        cat ~{write_lines(files)} | gsutil -m cp -I "~{gcs_output_dir}"
    >>>

    output {
        String gcs_dir = gcs_output_dir
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
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
