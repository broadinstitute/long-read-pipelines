version 1.0

##########################################################################################
# This WDL pipeline downloads files from wget-able URLs in parallel and stores the results in
# the specified GCS dir.  This pipeline is essentially a Cromwell/GCP reimagining of the
# Nextflow/AWS downloading pipeline from @alaincoletta (see: http://broad.io/aws_dl).
##########################################################################################

import "tasks/Structs.wdl"

workflow DownloadFromHA {
    input {
        String participant_id
        String sample_id
        String wget_rawdata_url
        String wget_md5_url

        String gcs_out_root_dir
    }

    parameter_meta {
        participant_id: "ID of the participant"
        sample_id:      "ID of the dataset belonging to the participant"
        wget_url:       "The wget-able URLs to download"
        gcs_output_dir: "GCS path for storing output"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call GetFileInfo { input: wget_url = wget_rawdata_url }

    # Assuming we have something in the manifest to download, download the files.
#    call DownloadFile {
#        input:
#            wget_url         = wget_url,
#            file_size_bytes  = GetFileInfo.file_size_bytes,
#            file_path        = GetFileInfo.file_path,
#            gcs_out_root_dir = outdir
#    }
}

# This task determines how much disk is needed for the part of the manifest that will be downloaded.
task GetFileInfo {
    input {
        String wget_url

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -x

        timeout --preserve-status 1 wget -S "~{wget_url}" 2> test.txt; grep '^Length' test.txt | awk '{ print $2 }' > size.txt
        echo '~{wget_url}' | sed 's/?.*//g' | awk -F"/" '{ print $NF }' > path.txt
    >>>

    output {
        Int file_size_bytes = read_int("size.txt")
        String file_path = read_string("path.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            1,
        boot_disk_gb:       10,
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

# This task checks if a file exists at the GCS output directory, and if not, downloads it from the server.
task DownloadFile {
    input {
        String wget_url
        Float file_size_bytes
        String file_path
        String gcs_out_root_dir

        RuntimeAttr? runtime_attr_override
    }

    # estimate required disk size in GB
    Float gb_bytes = 1024 * 1024 * 1024
    Int disk_size = 2 + 2*ceil((file_size_bytes / gb_bytes))

    command <<<
        set -x

        gcsfile="~{gcs_out_root_dir}/~{file_path}"

        if gsutil -q stat $gcsfile ; then
            echo "$gcsfile already exists."
        else
            wget -O ~{file_path} ~{wget_url}
            gsutil cp ~{file_path} $gcsfile
        fi
    >>>

    output {
        String out = read_string(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
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
