version 1.0

import "../../../structs/Structs.wdl"

workflow DownloadFromFTP {

    meta {
        description: "Download files from FTP in parallel and store the results in the specified GCS dir. This pipeline is essentially a Cromwell/GCP reimagining of the Nextflow/AWS downloading pipeline from @alaincoletta (see: http://broad.io/aws_dl)."
    }
    parameter_meta {
        ftp_dirs:                   "The FTP directories to download"
        num_simultaneous_downloads: "[default-valued] The number of files to fetch simultaneously."
        prepend_dir_name:           "If true, place the files in a subdirectory based on the basename of the FTP dir."
        exclude:                    "[default-valued] Simple substring patterns to exclude from the download."
        gcs_output_dir:             "GCS path for storing output"
    }

    input {
        Array[String] ftp_dirs

        Int num_simultaneous_downloads = 10
        Boolean prepend_dir_name = true
        Array[String] exclude = []

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    # Get a manifest of every single file to be downloaded, including file sizes.
    call GetFileManifest { input: ftp_dirs = ftp_dirs, exclude = exclude }

    # Parallelize downloads over num_simultaneous_downloads jobs.
    scatter (n in range(num_simultaneous_downloads)) {
        # Each of these jobs below will process (1/num_simultaneous_downloads)-th of the manifest.
        # We needn't actually shard the manifest beforehand; that'll be done within each of these
        # tasks using a nice awk idiom from Khalid Shakir.

        # Calculate the maximum amount of disk space the largest file in this manifest shard will need.
        call ComputeDiskSize {
            input:
                manifest = GetFileManifest.manifest,
                nth      = n,
                num_jobs = num_simultaneous_downloads,
        }

        # Assuming we have something in the manifest to download, download the files.
        if (ComputeDiskSize.max_size_bytes > 0) {
            call DownloadFTPFile {
                input:
                    manifest         = GetFileManifest.manifest,
                    nth              = n,
                    num_jobs         = num_simultaneous_downloads,
                    max_size_bytes   = ComputeDiskSize.max_size_bytes,
                    prepend_dir_name = prepend_dir_name,
                    gcs_out_root_dir = outdir
            }
        }
    }
}

# This task takes an FTP directory and emits a file list (path and sizes) to a table.
task GetFileManifest {
    input {
        Array[String] ftp_dirs
        Array[String] exclude

        RuntimeAttr? runtime_attr_override
    }

    Int num_excludes = length(exclude)

    command <<<
        set -euxo pipefail

        cat ~{write_lines(ftp_dirs)} > ftp_dirs.txt

        while read FTP_DIR
        do
            lftp -e 'find -l > all_files.txt; exit' $FTP_DIR
            grep -v '^d' all_files.txt | awk -v ftp_dir=$FTP_DIR '{ print ftp_dir, $6, $3 }' >> manifest.txt
        done < ftp_dirs.txt

        if [[ ~{num_excludes} -gt 0 ]]
        then
            grep -v -e ~{sep=' -e ' exclude} manifest.txt > manifest.new.txt
            mv manifest.new.txt manifest.txt
        fi
    >>>

    output {
        File manifest = "manifest.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            1,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-cloud-downloader:0.2.3"
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

# This task determines how much disk is needed for the part of the manifest that will be downloaded.
task ComputeDiskSize {
    input {
        File manifest
        Int nth
        Int num_jobs

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        echo 0 > sizes.txt
        awk 'NR % ~{num_jobs} == ~{nth} { print $3 }' ~{manifest} | while read ftp_size; do
            echo $ftp_size >> sizes.txt
        done

        cat sizes.txt

        sort -n sizes.txt | tail -1 > max_size.txt
    >>>

    output {
        Float max_size_bytes = read_float("max_size.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            1,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-cloud-downloader:0.2.3"
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

# This task checks if a file exists at the GCS output directory, and if not, downloads it from the FTP server.
task DownloadFTPFile {
    input {
        File manifest
        Int nth
        Int num_jobs
        Float max_size_bytes
        Boolean prepend_dir_name
        String gcs_out_root_dir

        RuntimeAttr? runtime_attr_override
    }

    # estimate required disk size in GB
    Float gb_bytes = 1024 * 1024 * 1024
    Int disk_size = 2 + ceil((max_size_bytes / gb_bytes))

    command <<<
        set -x

        RET=0

        awk 'NR % ~{num_jobs} == ~{nth} { print $0 }' ~{manifest} | while read line; do
            ftp_dir=$(echo $line | awk '{ print $1 }')
            ftp_base=$(basename $ftp_dir)
            ftp_file=$(echo $line | awk '{ print $2 }' | sed 's/^\.\///')

            gcsfile="~{gcs_out_root_dir}/~{true="$ftp_base/" false="" prepend_dir_name}$ftp_file"

            if gsutil -q stat "$gcsfile" ; then
                echo "$gcsfile already exists."
            else
                mkdir -p $(dirname $ftp_file)

                (lftp -e "get $ftp_file -o $ftp_file; exit" $ftp_dir && gsutil cp $ftp_file "$gcsfile") || RET=1

                rm -f $ftp_file
            fi
        done

        exit $RET
    >>>

    output {
        String out = read_string(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        3,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-cloud-downloader:0.2.3"
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
