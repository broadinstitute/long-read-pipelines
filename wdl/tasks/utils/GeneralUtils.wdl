## #

version 1.0

# todo: move all existing, simiar util tasks here
# Hosting utils that just needs basic Linux shell programs

task TarGZFiles {
    meta {
        description: "Zip up a list of files to a tar.gz file."
    }

    input {
        Array[File] files
        String name
    }

    command <<<
        set -eux
        mkdir -p save/
        for ff in ~{sep=' ' files}; do cp "${ff}" save/ ; done
        tar -cvzf ~{name}.tar.gz -C save/ .
    >>>

    output {
        File you_got_it = "~{name}.tar.gz"
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task GetTodayDate {
    meta {
        description: "Generates a YYYY-MM-DD date of today (when this task is called). UTC."
        volatile: true
    }
    command {
        date '+%Y-%m-%d'
    }

    output {
        String yyyy_mm_dd = read_string(stdout())
    }
    runtime {
        disks: "local-disk 10 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}
