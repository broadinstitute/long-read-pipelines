version 1.0

workflow TestBatchAttributes {

    meta {
        author: "Jonn Smith"
        description: "This workflow tests the batch attributes"
    }

    call GetCurrentTimestampString as t_001_WdlExecutionStartTimestamp { input: }

    output {
        String timestamp_string = t_001_WdlExecutionStartTimestamp.timestamp_string
    }

}


task GetCurrentTimestampString {

    meta {
        # The volatile keyword forces call caching to be turned off, which is
        # exactly what we want for this task.
        # For more info see: https://cromwell.readthedocs.io/en/stable/optimizations/VolatileTasks/
        volatile: true
        description: "Get the current timestamp as a string"
    }

    parameter_meta {
        date_format: "The date format string to use. See the unix `date` command for more info."
    }

    input {
        String date_format = "%Y%m%d_%H%M%S_%N"
    }

    String date_file = "the_date_file.txt"

    command {
        date +~{date_format} > ~{date_file}
        cat ~{date_file}
    }

    # ------------------------------------------------
    # Runtime settings:
     runtime {
         docker: "ubuntu:19.10"
         memory: "512 MB"
         disks: "local-disk 10 HDD"
         bootDiskSizeGb: 15
         preemptible: 0
         cpu: 1
     }

    output {
        String timestamp_string   = read_string(date_file)
    }
}
