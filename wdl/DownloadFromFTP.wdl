version 1.0

# Copyright Broad Institute, 2019
#
# About:
#   This WDL pipeline downloads data from FTPs in parallel and stores the results in the
#   specified GCS dir.  This pipeline is essentially a Cromwell/GCP reimagining of the
#   Nextflow/AWS downloading pipeline from @alaincoletta (see: http://broad.io/aws_dl).
#
# Description of inputs:
#   Required:
#       Array[String] URLS                    - The URL(s) of the datasets to download.
#       String gcs_output_dir                 - GCS output dir.
#
# Licensing:
#   This script is released under the WDL source code license (BSD-3) (see LICENSE in
#   https://github.com/broadinstitute/wdl). Note however that the programs it calls may
#   be subject to different licenses. Users are responsible for checking that they are
#   authorized to run all programs before running this script.

workflow DownloadFromFTP {
    input {
        Array[String] URLs
        String gcs_output_dir
    }

    scatter (URL in URLs) {
        call Download {
            input:
                URL            = URL,
                gcs_output_dir = gcs_output_dir
        }
    }
}

# This task checks to see if a to-be-downloaded file exists at the specified GCS filepath,
# and if not, initiates a download process.  The downloaded file is then immediately
# uploaded to the specified GCS filepath.  It is then # deleted from local storage here by
# Cromwell, preventing redundant storage of the data.

task Download {
    input {
        String URL
        String gcs_output_dir
    }

    Int cpus = 2
    String basename = basename(URL)
    String out_path = sub(gcs_output_dir, "/$", "") + "/" + basename

    command <<<
        set -euxo pipefail

        if gsutil ls -lh ~{out_path} ; then
            echo "~{out_path} already exists."
        else
            wget -q ~{URL}
            gsutil -o "GSUtil:parallel_process_count=~{cpus}" -o "GSUtil:parallel_thread_count=1" -m cp ~{basename} ~{gcs_output_dir}
        fi

        gsutil ls -lh ~{out_path} > upload_list.txt
    >>>

    output {
        File upload_list = "upload_list.txt"
    }

    runtime {
        cpu: "~{cpus}"
        memory: "4GB"
        disks: "local-disk 300 LOCAL"
        preemptible: 0
        maxRetries: 0
        docker: "kgarimella/lr-cloud-downloader:0.02.0"
    }
}
