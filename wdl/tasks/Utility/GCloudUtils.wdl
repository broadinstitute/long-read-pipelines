version 1.0

task GenerateDummyFile {
    meta {
        desciption:
        "Generate a dummy file. Currently the only purpose of this task is to assist a trick that allows manual file delocalization in a task."
    }
    output {
        File dummy = "dummy.txt"
    }

    command <<<
    set -euxo pipefail
        echo "don't touch me" > "dummy.txt"
    >>>
    runtime {
        preemptible_tries:     1
        max_retries:           1
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task GetParentDir {
    meta {
        desciption:
        "Get the parent directory of a given GCS path"
    }
    input {
        String path
    }

    output {
        String parent_path = read_string("results.txt")
    }

    String pp = sub(path, "/+$", "")

    command <<<
    set -euxo pipefail
        parent_path=$(dirname ${pp})
        echo "${parent_path}" > results.txt
    >>>
    runtime {
        preemptible_tries:     1
        max_retries:           1
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task ListFilesInDir {
    meta {
        desciption:
        "List files in a given GCS 'directory'"
    }
    input {
        String path
        String? pattern
    }

    output {
        Array[String] files = read_lines("results.txt")
    }

    String compile_time_pattern = select_first([pattern, ".*"])

    command <<<
    set -euxo pipefail
        if ~{defined(pattern)}; then
            gsutil ls ${path} | grep -E "~{compile_time_pattern}" > results.txt
        else
            gsutil ls ${path} > results.txt
        fi
    >>>

    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 10 HDD"
        preemptible_tries:     1
        max_retries:           1
        docker:"us.gcr.io/google.com/cloudsdktool/google-cloud-cli:alpine"
    }
}