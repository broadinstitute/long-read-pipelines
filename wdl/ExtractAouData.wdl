version 1.0

##########################################################################################
# This WDL pipeline downloads files from wget-able URLs in parallel and stores the results in
# the specified GCS dir.  This pipeline is essentially a Cromwell/GCP reimagining of the
# Nextflow/AWS downloading pipeline from @alaincoletta (see: http://broad.io/aws_dl).
##########################################################################################

import "tasks/Structs.wdl"

workflow ExtractAouData {
    input {
        File gs_path
        File gs_md5

        String gcs_out_root_dir
    }

    parameter_meta {
        gs_path:          "Path to tarball"
        gs_md5:           "Path to MD5 file for tarball"
        gcs_out_root_dir: "GCS path for storing output"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call VerifyAndExtractTarball {
        input:
            gs_path          = gs_path,
            gs_md5           = gs_md5,
            gcs_out_root_dir = outdir
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
    String bn = sub(basename(gs_path), "_[1234]_[ABCD]0[1234]_rawdata.tar.gz", "")

    command <<<
        set -euxo pipefail

        echo ~{disk_size}
        echo ~{bn}
        du -hcs .

#        mkdir -p /dls/storage/longreadtemp/
#        mv ~{gs_path} /dls/storage/longreadtemp/
#        md5sum --check ~{gs_md5}
#
#        tar zxvf ~{gs_path}
#
#        find . -exec ls -lah {} \;
#        find . \( -name \*.bam -or -name \*.pbi -or -name \*.xml \) mv {} \.
#
#        ls -lah
    >>>

    output {
        File bam = "~{bn}.subreads.bam"
        File pbi = "~{bn}.subreads.bam.pbi"
        File xml = "~{bn}.subreadsset.xml"
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
