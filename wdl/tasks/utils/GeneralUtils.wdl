version 1.0

task GetTodayDate {
    meta {
        desciption: "Generates a YYYY-MM-DD date of today (when this task is called). UTC."
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
