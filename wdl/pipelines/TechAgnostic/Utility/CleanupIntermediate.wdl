version 1.0

workflow CleanupIntermediate {
    # Ironicaly, this generates intermeidate files too, but they are tiny.
    meta {
        description: "A workflow to clean up intermediate files from running workflows on Terra. Use at your own risk."
    }
    parameter_meta {
        submissionIDs: "List of submissions whose intermediate files are to be deleted"
    }

    input {
        String workspace_bucket
        Array[String] submissionIDs
        Boolean keep_logs = false
    }

    scatter (sid in submissionIDs) {
        call CleanupAFolder { input:
            bucket_name = workspace_bucket,
            submission_id = sid,
            keep_logs = keep_logs
        }
    }
}

task CleanupAFolder {
    input {
        String bucket_name
        String submission_id
        Boolean keep_logs = false
    }

    Boolean fail = submission_id == ""
    command <<<
        if ~{fail}; then exit 1; fi
        if ~{keep_logs}; then
            # keep some (presumably) lightweight cromwell log files
            timeout 23h \
            gcloud -q storage ls --recursive \
                gs://~{bucket_name}/submissions/~{submission_id} \
                | grep -v "*log$" | grep -vF '/stderr' | grep -vF '/stdout' | grep -vF '/script' \
                | gcloud -q storage rm -rf -I \
            || echo "Timed out. Please try again."
        else
            timeout 23h \
            gcloud -q storage rm \
                -rf gs://~{bucket_name}/submissions/~{submission_id} \
            || echo "Timed out. Please try again."
        fi
    >>>

    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 10 HDD"
        preemptible: 5
        maxRetries:  1
        docker:"us.gcr.io/google.com/cloudsdktool/google-cloud-cli:alpine"
    }
}
