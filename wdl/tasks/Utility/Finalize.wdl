version 1.0

import "../../structs/Structs.wdl"

task FinalizeToFile {

    meta{
        description: "Copies the given file to the specified bucket."
    }

    parameter_meta {
        file: {
            description: "file to finalize",
            localization_optional: true
        }
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
        outdir: "directory to which files should be uploaded"
        name:   "name to set for uploaded file"
    }

    input {
        File file
        String outdir
        String? name

        File? keyfile

        RuntimeAttr? runtime_attr_override
    }



    String gcs_output_dir = sub(outdir, "/+$", "")
    String gcs_output_file = gcs_output_dir + "/" + select_first([name, basename(file)])

    command <<<
        set -euxo pipefail

        gsutil -m cp "~{file}" "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       25,
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

task FinalizeToDir {

    meta {
        description: "Copies the given file to the specified bucket."
    }

    parameter_meta {
        files: {
            description: "files to finalize",
            localization_optional: true
        }
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
        outdir: "directory to which files should be uploaded"
    }

    input {
        Array[File] files
        String outdir

        File? keyfile

        RuntimeAttr? runtime_attr_override
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
        boot_disk_gb:       25,
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

task FinalizeTarGzContents {
    meta {
        description : "Copies the contents of the given tar.gz file to the specified bucket."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        tar_gz_file : "Gzipped tar file whose contents we'll copy."
        outdir : "Google cloud path to the destination folder."

        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."

        runtime_attr_override : "[optional] Additional runtime parameters."
    }

    input {
        File tar_gz_file
        String outdir

        File? keyfile

        RuntimeAttr? runtime_attr_override
    }

    # This idiom ensures that we don't accidentally have double-slashes in our GCS paths
    String gcs_output_dir = sub(sub(outdir + "/", "/+", "/"), "gs:/", "gs://")

    command <<<
        set -euxo pipefail

        mkdir tmp
        cd tmp
        tar -zxf ~{tar_gz_file}

        gsutil -m cp -r * ~{gcs_output_dir}
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            10,
        boot_disk_gb:       25,
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

task WriteCompletionFile {

    meta {
        description : "Write a file to the given directory indicating the run has completed."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        outdir : "Google cloud path to the destination folder."
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
    }

    input {
        String outdir
        File? keyfile
    }

    command <<<
        set -euxo pipefail

        completion_file="COMPLETED_AT_$(date +%Y%m%dT%H%M%S).txt"
        touch $completion_file

        gsutil cp $completion_file ~{outdir}
    >>>

    #########################

    runtime {
        cpu:                    1
        memory:                 1 + " GiB"
        disks: "local-disk " +  10 + " HDD"
        bootDiskSizeGb:         25
        preemptible:            2
        maxRetries:             2
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
}

task CompressAndFinalize {

    meta {
        description: "Gzip a file and finalize"
    }

    parameter_meta {
        file : "File to compress and finalize."
        outdir : "Google cloud path to the destination folder."
        name : "[optional] Name of the file to write.  If not specified, the name of the input file will be used."
        runtime_attr_override : "[optional] Additional runtime parameters."
    }

    input {
        File file
        String outdir
        String? name

        RuntimeAttr? runtime_attr_override
    }

    String base = basename(file)
    String out = sub(select_first([name, base]), ".gz$", "") +  ".gz"
    # THIS IS ABSOLUTELY CRITICAL: DON'T CHANGE TYPE TO FILE, AS CROMWELL WILL TRY TO LOCALIZE THIS NON-EXISTENT FILE
    String gcs_output_file = sub(outdir, "/+$", "") + "/" + out

    Int disk_size = 2 * ceil(size(file, "GB"))

    command <<<
        set -euxo pipefail

        gzip -vkc ~{file} > "~{base}.gz"
        gsutil cp "~{base}.gz" "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
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

task FinalizeAndCompress {
    meta {
        description: "Gzip a bunch of files and finalize to the same \'folder\'"
    }

    parameter_meta {
        files : "Files to compress and finalize."
        outdir : "Google cloud path to the destination folder."
        prefix : "[optional] Prefix to add to the output files."
        runtime_attr_override : "[optional] Additional runtime parameters."
    }

    input {
        Array[File] files
        String outdir

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    String gcs_output_file = sub(outdir, "/+$", "") + "/" + prefix + "/"

    Int disk_size = 5 * ceil(size(files, "GB"))

    command <<<
        set -euxo pipefail

        for ff in ~{sep=' ' files};
        do
            base="$(basename -- ${ff})"
            mv "${ff}" "${base}" && gzip -vk "${base}"
        done

        gsutil -m cp /cromwell_root/*.gz "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             7,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
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

task WriteNamedFile {

    meta {
        description : "Write a file to the given directory with the given name."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String name
        String outdir
        File? keyfile
    }

    parameter_meta {
        name : "Name of the file to write."
        outdir : "Google cloud path to the destination folder."
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
    }

    command <<<
        set -euxo pipefail

        touch "~{name}"

        gsutil cp "~{name}" ~{outdir}
    >>>

    #########################

    runtime {
        cpu:                    1
        memory:                 1 + " GiB"
        disks: "local-disk " +  10 + " HDD"
        bootDiskSizeGb:         25
        preemptible:            2
        maxRetries:             2
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
}
