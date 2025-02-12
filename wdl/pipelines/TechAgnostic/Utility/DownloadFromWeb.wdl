version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils

workflow DownloadFromWeb {

    meta {
        description: "This WDL pipeline downloads directories from HTTP/FTP/SFTP servers in parallel and stores the results in the specified GCS dir.  This pipeline is essentially a Cromwell/GCP reimagining of the Nextflow/AWS downloading pipeline from @alaincoletta (see: http://broad.io/aws_dl)."
    }

    parameter_meta {
        manifest:                   "A file with a list of SRA ID(s) to download on each line"
        num_simultaneous_downloads: "[default-valued] The number of files to fetch simultaneously."
        prepend_dir_name:           "If true, place the files in a subdirectory based on the basename of the FTP dir."
        gcs_out_root_dir:           "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        File manifest

        Int num_simultaneous_downloads = 10
        Boolean prepend_dir_name = true

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    # Parallelize downloads over num_simultaneous_downloads jobs.
    scatter (n in range(num_simultaneous_downloads)) {
        # Each of these jobs below will process (1/num_simultaneous_downloads)-th of the manifest.
        # We needn't actually shard the manifest beforehand; that'll be done within each of these
        # tasks using a nice awk idiom from Khalid Shakir.

        # Assuming we have something in the manifest to download, download the files.
        call DownloadFiles {
            input:
                manifest         = manifest,
                nth              = n,
                num_jobs         = num_simultaneous_downloads,
                prepend_dir_name = prepend_dir_name,
                gcs_out_root_dir = outdir
        }
    }
}

task DownloadFiles {

    meta {
        description: "This task checks to see if a to-be-downloaded file exists at the specified GCS filepath, and if not, initiates a parallel download process.  The downloaded file is then immediately uploaded to the specified GCS filepath.  It is then deleted from local storage here by Cromwell, preventing redundant storage of the data."
    }

    input {
        File manifest
        Int nth
        Int num_jobs
        Boolean prepend_dir_name
        String gcs_out_root_dir

        Int disk_size_gb = 100
        Int num_cpus = 4

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -x

        RET=0

        awk 'NR % ~{num_jobs} == ~{nth} { print $0 }' ~{manifest} | while read line; do
            name=$(echo $line | xargs -n 1 basename)
            dir=$(echo $line | perl -pe 's|(^.*://.*?/)||' | xargs -n 1 dirname)
            gcsdir="~{gcs_out_root_dir}~{true="/$dir" false="" prepend_dir_name}"

            if gsutil -q stat "${gcsdir}/${name}" ; then
                echo "${gcsdir}/${name} already exists."
            else
                (aria2c $line && gsutil cp $name "$gcsdir/") || RET=1

                rm -f $name
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
        disk_gb:            disk_size_gb,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        3,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-cloud-downloader:0.2.5"
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
