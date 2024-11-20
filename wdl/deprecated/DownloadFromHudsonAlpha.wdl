version 1.0

##########################################################################################
# This WDL pipeline downloads files from wget-able URLs in parallel and stores the results in
# the specified GCS dir.  This pipeline is essentially a Cromwell/GCP reimagining of the
# Nextflow/AWS downloading pipeline from @alaincoletta (see: http://broad.io/aws_dl).
##########################################################################################

import "../structs/Structs.wdl"

workflow DownloadFromHudsonAlpha {
    input {
        File manifest

        String gcs_out_root_dir
    }

    parameter_meta {
        urls:             "Files to download"
        gcs_out_root_dir: "GCS path for storing output"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    scatter (l in read_tsv(manifest)) {
        call GetFileSize {
            input:
                url = l[0],
                temp_url_sig = l[1],
                temp_url_expires = l[2]
        }

        call DownloadFile {
            input:
                url = l[0],
                temp_url_sig = l[1],
                temp_url_expires = l[2],
                size = GetFileSize.size
        }

        call ExtractFiles {
            input:
                file = DownloadFile.outfile,
                gcs_out_root_dir = outdir
        }
    }
}

task GetFileSize {
    input {
        String url
        String temp_url_sig
        String temp_url_expires

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -x

        timeout 2 wget -O /dev/null "~{url}?temp_url_sig=~{temp_url_sig}&temp_url_expires=~{temp_url_expires}" 2> header.txt
        grep 'Length' header.txt | awk '{ print int(($2/1e9)+1) }' > size.txt
    >>>

    output {
        Int size = read_int("size.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            1,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        2,
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

task DownloadFile {
    input {
        String url
        String temp_url_sig
        String temp_url_expires
        Int size

        RuntimeAttr? runtime_attr_override
    }

    String fn = basename(url)

    command <<<
        set -euxo pipefail

        wget -O "~{fn}" "~{url}?temp_url_sig=~{temp_url_sig}&temp_url_expires=~{temp_url_expires}"
    >>>

    output {
        File outfile = fn
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            size,
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

task VerifyAndExtractTarball {
    input {
        File gs_path
        File gs_md5
        String gcs_out_root_dir

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size([gs_path, gs_md5], "GB"))
    String bn = sub(basename(gs_path), "_rawdata.tar.gz", "")

    command <<<
        set -euxo pipefail

        tar zxvf ~{gs_path}

        find . -exec ls -lah {} \;

        find . \( -name \*.bam -or -name \*.pbi -or -name \*.xml \) \
            -exec gsutil cp {} "~{gcs_out_root_dir}/inputs/~{bn}/" \;

        gsutil ls "~{gcs_out_root_dir}/inputs/~{bn}/*.bam" > bam.txt
        gsutil ls "~{gcs_out_root_dir}/inputs/~{bn}/*.pbi" > pbi.txt
        gsutil ls "~{gcs_out_root_dir}/inputs/~{bn}/*.xml" > xml.txt
    >>>

    output {
        String bam = read_string("bam.txt")
        String pbi = read_string("pbi.txt")
        String xml = read_string("xml.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
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

task ExtractFiles {
    input {
        File file
        String gcs_out_root_dir

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(file, "GB"))

    command <<<
        set -euxo

        mkdir out

        if [[ "~{file}" =~ \.tar.gz$ ]]; then
            tar zxvf ~{file} -C out/ || tar xvf ~{file} -C out/
        else
            mv ~{file} out/~{basename(file)}
        fi

        cd out
        gsutil -m cp -r * ~{gcs_out_root_dir}/

        cd ..
        touch done
    >>>

    output {
        File done = "done"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
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
