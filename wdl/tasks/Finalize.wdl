version 1.0

import "Structs.wdl"

task FinalizeToDir {
    input {
        Array[String] files
        String outdir
    }

    # This idiom ensures that we don't accidentally have double-slashes in our GCS paths
    String gcs_output_dir = sub(sub(outdir + "/", "/+", "/"), "gs:/", "gs://")

    command <<<
        set -euxo pipefail

        gsutil -m cp ~{sep=' ' files} ~{gcs_output_dir}
    >>>

    output { }

    runtime {
        cpu:                    1
        memory:                 2
        disks:                  "local-disk 10 HDD"
        bootDiskSizeGb:         10
        preemptible:            2
        maxRetries:             2
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
}

