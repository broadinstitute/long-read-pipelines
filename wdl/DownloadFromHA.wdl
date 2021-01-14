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

    call DownloadFile {
        input:
            wget_rawdata_url = wget_rawdata_url,
            gcs_out_root_dir = outdir
    }
}

task DownloadFile {
    input {
        String wget_rawdata_url
        String gcs_out_root_dir

        RuntimeAttr? runtime_attr_override
    }

    String file_path = basename(sub(wget_rawdata_url, "\?.*", ""))
    Int disk_size = 1

    command <<<
        set -euxo pipefail

        #wget -O /mnt/data/~{file_path} ~{wget_rawdata_url}
        echo ~{file_path}
        echo ~{wget_rawdata_url}
        ls /mnt/data/
    >>>

    output {
        File out = stdout()
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-cloud-downloader:0.2.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "/mnt/data " +   select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
