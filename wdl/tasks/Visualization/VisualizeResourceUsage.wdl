version 1.0

task SimpleRscript {
    meta {
        description: "Uses a simple R script to visualize resource usage of a resource-hungry task."
    }
    parameter_meta {
        resource_log: "Resource usage log file."
        output_pdf_name: "name of the output plot, must end in .pdf"
        plot_title: "Title of the plot"
        format: "The format that the monitor log represent time in, e.g. 'Sat Oct 28 01:49:01 UTC 2023' would use '%a %b %d %H:%M:%S %Y', and 'Fri 10 Nov 2023 01:31:12 AM UTC' would use '%a %d %b %Y %H:%M:%S'"
        use_us_eastern_tz: "if true, display time in the plot using the Boston timezone"
    }
    input {
        File resource_log
        String output_pdf_name
        String plot_title

        String format = '%Y-%m-%d-%H-%M-%S %a'
        Boolean use_us_eastern_tz = false
    }
    output {
        File plot_pdf = "~{output_pdf_name}"
    }

    command <<<
        set -eux

        /opt/plot.resources.R \
            -I "~{resource_log}" \
            -O "~{output_pdf_name}" \
            --title "~{plot_title}" \
            --format "~{format}" \
            ~{true='--boston' false=' ' use_us_eastern_tz}

    >>>
    runtime {
        disks: "local-disk 10 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-resource-visual:0.1.2"
    }
}
