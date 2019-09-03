version 1.0

# TODO: describe purpose
task DetectRunInfo {
    input {
        String gcs_dir
        String? sample_name
        String docker
    }

    String SM = if defined(sample_name) then "--SM " + sample_name else ""

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        python /usr/local/bin/detect_run_info.py ~{SM} ~{gcs_dir} > run_info.txt
        gsutil ls ~{gcs_dir} | grep -v scraps | grep -e '\.bam$' -e '\.f\(ast\)\?q\(\.gz\)\?$' > files.txt
    >>>

    output {
        File run_info_file = "run_info.txt"
        File fofn = "files.txt"
        Map[String, String] run_info = read_map("run_info.txt")
        Array[String] files = read_lines("files.txt")
    }

    runtime {
        cpu: 1
        memory: "1GB"
        disks: "local-disk 1 SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}

# TODO: describe purpose
task PrepareRun {
    input {
        Array[File] files
        String docker
    }

    Int disk_size = 2*ceil(size(files, "GB"))

    command <<<
        set -euxo pipefail

        python /usr/local/bin/prepare_run.py ~{sep=' ' files}
    >>>

    output {
        File unmapped_bam = "unmapped.bam"
    }

    runtime {
        cpu: 2
        memory: "2GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}