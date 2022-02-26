version 1.0

workflow CleanupIntermediate {
    # Ironicaly, this generates intermeidate files too, but they are tiny.
    meta {
        description: "A workflow to clean up intermediate files from running workflows. Use at your own risk."
    }

    input {
        String workspace_bucket
        Array[String] submissionIDs
    }

    parameter_meta {
        submissionIDs: "List of submissions whose intermediate files are to be deleted"
    }

    scatter (sid in submissionIDs) {
        call CleanupAFolder {
            input:
                bucket_name = workspace_bucket,
                submission_id = sid
        }
    }
}

task CleanupAFolder {
    input {
        String bucket_name
        String submission_id
    }

    command <<<
        echo "started"
        gsutil -q rm -rf gs://~{bucket_name}/~{submission_id}
    >>>

    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           1
        docker:"google/cloud-sdk:latest"
    }
}