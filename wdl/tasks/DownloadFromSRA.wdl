version 1.0

# Copyright Broad Institute, 2019
#
# About:
#   This WDL pipeline downloads data from SRA in parallel and stores the results in the
#   specified GCS dir.  This pipeline is essentially a Cromwell/GCP reimagining of the
#   Nextflow/AWS downloading pipeline from @alaincoletta (see: http://broad.io/aws_dl).
#
# Description of inputs:
#   Required:
#       Array[String] SRA_IDs                 - The SRA ID(s) of the datasets to download.
#       String gcs_output_dir                 - GCS output dir.
#
# Licensing:
#   This script is released under the WDL source code license (BSD-3) (see LICENSE in
#   https://github.com/broadinstitute/wdl). Note however that the programs it calls may
#   be subject to different licenses. Users are responsible for checking that they are
#   authorized to run all programs before running this script.

workflow DownloadFromSRA {
    input {
        Array[String] SRA_IDs
        String gcs_output_dir
    }

    scatter (SRA_ID in SRA_IDs) {
        call GetSRAIDs {
            input:
                SRA_ID                         = SRA_ID
        }

        scatter (row in read_tsv(GetSRAIDs.sra_table)) {
            call FastqDump {
                input:
                    SRR_ID                     = row[0],
                    gcs_output_dir             = gcs_output_dir,
                    uncompressed_disk_space_mb = row[1]
            }
        }
    }
}

# This task takes an SRA ID of some sort (e.g. study ID, run ID, experiment ID, etc.), finds the SRR/ERR run IDs
# associated with the provided SRA ID, and writes their accession numbers and file sizes to a tsv table.
#
# In inspecting the filesizes returned (column 8), it appears empirically to be around 5 times larger than the
# gzipped fastq files themselves.  Thus I assume this number is the uncompressed size of the file, but I cannot
# find any documentation that explicitly states this.  I've assumed a lesser compression level of 2-fold below.

task GetSRAIDs {
    input {
        String SRA_ID
    }

    command <<<
        esearch -db sra -query ~{SRA_ID} | efetch --format runinfo | grep '[SE]RR' | cut -f 1,8 -d',' | sed 's/,/\t/g' > sra.txt
    >>>

    output {
        File sra_table = "sra.txt"
    }

    runtime {
        cpu: 1
        memory: "2GB"
        disks: "local-disk 1 LOCAL"
        preemptible: 1
        maxRetries: 0
        docker: "us.gcr.io/broad-dsp-lrma/lr-cloud-downloader:0.2.1"
    }
}

# This task checks to see if a to-be-downloaded file exists at the specified GCS filepath, and if not, initiates a
# parallel download process.  The downloaded file is then immediately uploaded to the specified GCS filepath.  It
# is then deleted from local storage here by Cromwell, preventing redundant storage of the data.
#
# In the future, it might be interesting to explore gcsfuse (https://cloud.google.com/storage/docs/gcs-fuse) as an
# option for downloading data directly to a GCS filepath.
task FastqDump {
    input {
        String SRR_ID
        String gcs_output_dir
        Int uncompressed_disk_space_mb
    }

    Int cpus = 4

    Int mb_to_gb_divisor = 1000
    Int compression_factor = 2
    Int disk_space_divisor = mb_to_gb_divisor * compression_factor
    Int minimum_disk_size_gb = 10
    Int compressed_disk_space_gb = ceil(uncompressed_disk_space_mb / disk_space_divisor)

    # We're going to request a LOCAL disk which comes in minimum increments of 375 GB and ignores our disk space
    # request.  But we're going to compute it anyway should we decide to change our request to an HDD/SSD later.
    Int final_disk_space_gb = if (compressed_disk_space_gb < minimum_disk_size_gb) then minimum_disk_size_gb else compressed_disk_space_gb

    command <<<
        if gsutil ls -lh ~{gcs_output_dir}/~{SRR_ID}.fastq.gz ; then
            echo "~{gcs_output_dir}/~{SRR_ID} already exists."
        else
            set -eu # we place this here because the gsutil check up above will error us out were we to place it at the very beginning

            # we do these because sratools actually prefetch SRA files to $HOME,
            # which in cromwell-controlled workflows typically are quite small, and you'll see errors like
            # "fasterq-dump.2.9.1 int: storage exhausted while writing file within file system module"
            # then ultimately may fail to out of space errors
            mkdir -p $HOME/.ncbi/ && mkdir -p ncbi/public/ && \
                echo '/repository/user/main/public/root = "/cromwell_root/ncbi/public"' > $HOME/.ncbi/user-settings.mkfg

            echo "=========="; df -h; echo "==========";
            fasterq-dump ~{SRR_ID} --threads ~{cpus}
            du -sh *
            echo "-----"
            for fastq in *.fastq; do pigz $fastq; done
            echo "-----"
            du -sh *
            echo "=========="; df -h; echo "==========";
            gsutil -o "GSUtil:parallel_process_count=~{cpus}" \
                   -o "GSUtil:parallel_thread_count=1" \
                   -m cp *.gz ~{gcs_output_dir}
            gsutil ls -lh ~{gcs_output_dir} > upload_list.txt
        fi

        exit 1
    >>>

    output {
        File upload_list = "upload_list.txt"
    }

    runtime {
        cpu: cpus
        memory: "4GB"
        disks: "local-disk ~{final_disk_space_gb} LOCAL"
        preemptible: 0
        maxRetries: 0
        docker: "us.gcr.io/broad-dsp-lrma/lr-cloud-downloader:0.2.1"
    }
}
