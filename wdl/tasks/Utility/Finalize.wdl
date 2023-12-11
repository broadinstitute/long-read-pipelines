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

        gcloud storage cp "~{file}" "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
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
        file_names: "custom names for files; must be the same length as files if provided"
        outdir: "directory to which files should be uploaded"

        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
    }

    input {
        Array[File] files
        Array[String]? file_names
        String outdir

        File? keyfile

        RuntimeAttr? runtime_attr_override
    }

    String gcs_output_dir = sub(outdir, "/+$", "")

    Boolean fail = if(defined(file_names)) then length(select_first([file_names])) == length(files) else false
    # this variable is defined because of meta-programing:
    # Cromwell generates the script to be executed at runtime (duing the run of the workflow),
    # but also at "compile time" when looked from the individual task perspective--the task is "compiled" right before it is run.
    # so optional variables, if not specified, cannot be used in the command section because at that "compile time", they are undefined
    # here we employ a hack:
    # if the optional input file_names isn't provided, it's not used anyway, so we don't worry about the literal correctness of
    # the variable's values--the variable used in generating the script--but only care that it is defined.
    Array[String] names_for_cromwell = select_first([file_names, ["correctness_doesnot_matter_here"]])
    command <<<
        set -euxo pipefail

        if ~{fail}; then echo "input files and file_names don't have the same length!" && exit 1; fi

        if ~{defined(file_names)}; then
            paste \
                ~{write_lines(files)} \
                ~{write_lines(names_for_cromwell)} \
            > file_and_customname.tsv
            while IFS=$'\t' read -r ff nn; do
                gcloud storage cp \
                    "${ff}" \
                    "~{gcs_output_dir}"/"${nn}"
            done < file_and_customname.tsv
        else
            cat ~{write_lines(files)} | \
            gcloud storage cp -I "~{gcs_output_dir}"
        fi
    >>>

    output {
        String gcs_dir = gcs_output_dir
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
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
        bootDiskSizeGb:         10
        preemptible:            2
        maxRetries:             2
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
}

task WriteNamedFile {

    meta {
        description : "Write a file to the given directory with the given name."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        name : "Name of the file to write."
        outdir : "Google cloud path to the destination folder."
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
    }

    input {
        String name
        String outdir
        File? keyfile
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
        bootDiskSizeGb:         10
        preemptible:            2
        maxRetries:             2
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
}

task CompressAndFinalize {

    meta {
        description: "(Block-)Gzip a file and finalize"
    }

    parameter_meta {
        file : {desciption: "File to compress and finalize.", localization_optional: true}
        block_gzip: "if true, will block-gzip the file (preferrable for certain genomics files)"

        outdir : "Google cloud path to the destination folder."
        name : "[optional] Name of the file to write.  If not specified, the name of the input file will be used."
        runtime_attr_override : "[optional] Additional runtime parameters."
    }

    input {
        File file
        String outdir
        String? name

        Boolean block_gzip = false
        String disk_type = "SSD"
        RuntimeAttr? runtime_attr_override
    }

    String base = basename(file)
    String out = sub(select_first([name, base]), ".gz$", "") +  ".gz"

    # THIS IS ABSOLUTELY CRITICAL: DON'T CHANGE TYPE TO FILE, AS CROMWELL WILL TRY TO LOCALIZE THIS NON-EXISTENT FILE
    String gcs_output_file = sub(outdir, "/+$", "") + "/" + out

    command <<<
        set -euxo pipefail

        time \
        gcloud storage cp ~{file} localized

        if ~{block_gzip}; then
            bgzip -kc -t ~{cores} localized > "~{base}.gz"
        else
            pigz -vkc localized > "~{base}.gz"
        fi
        gsutil cp "~{base}.gz" "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    Int disk_size = 2 * ceil(size(file, "GiB"))
    Int cores = if (disk_size>4) then 4 else 1
    RuntimeAttr default_attr = object {
        cpu_cores:          cores,
        mem_gb:             2*cores,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
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
        files : {desciption: "Files to compress and finalize.", localization_optional: true}
        outdir : "Google cloud path to the destination folder."
        folder : "new folder under 'outdir' to hold the compressed files."
        runtime_attr_override : "[optional] Additional runtime parameters."
    }

    input {
        Array[File] files
        Boolean block_gzip = false

        String outdir
        String? folder

        String disk_type = "HDD"
        RuntimeAttr? runtime_attr_override
    }

    String gcs_output_dir = sub(outdir, "/+$", "") + "/" + if (defined(folder)) then "~{folder}/" else ""

    command <<<
        set -euxo pipefail

        mkdir -p copy
        time \
        gcloud storage cp ~{sep=' ' files} copy/

        cd copy
        for ff in *; do
            if ~{block_gzip}; then bgzip -k "${ff}" ; else pigz -vk "${ff}"; fi
        done
        ls &&
        cd - &&
        gcloud storage rsync \
            copy/ \
            "~{gcs_output_dir}"
    >>>

    output {
        String gcs_path = gcs_output_dir
    }

    #########################
    Int disk_size = 2 * ceil(size(files, "GiB"))

    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task TarGZFilesAndSave {
    meta {
        desciption:
        ""
    }

    parameter_meta {
        files: {localization_optional: true}
        name: ""
    }

    input {
        Array[File]+ files
        String name
        String outdir

        String disk_type = "SSD"
        RuntimeAttr? runtime_attr_override
    }

    output {
        String gcs_path = gcs_output_file
    }

    String gcs_output_dir = sub(outdir, "/+$", "")
    String local_out = "~{name}"
    String gcs_output_file = gcs_output_dir + "/~{local_out}"

    command <<<
        set -euxo pipefail

        mkdir -p local
        time \
        gcloud storage cp ~{sep=' ' files} \
        local/

        tar cf - local/* | \
        pigz \
        > "~{local_out}"

        time \
        gcloud storage cp "~{local_out}" "~{gcs_output_file}"
    >>>

    #########################
    Int disk_size = 2 * ceil(size(files, "GiB"))

    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
