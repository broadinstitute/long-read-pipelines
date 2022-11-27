version 1.0

task GetBlobFolder {
    meta {
        description: "Getting the bucket and folder of a GCS object (eg. gs://bucket/folder/my_object.txt -> gs://bucket/folder)"
    }

    input {
        String blob
    }
    parameter_meta {
        blob: "A GCS path representing a blob/object."
    }

    command <<<
        set -eux

        echo ~{blob} | awk -F '/' 'BEGIN{OFS="/"} NF{NF-=1};1' > prefix.txt
    >>>

    output {
        String bucket_and_folder = read_string("prefix.txt")
    }

    runtime {
        cpu: 1
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
