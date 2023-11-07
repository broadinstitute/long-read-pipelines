version 1.0

task SimpleRscript {
    meta {
        description: "Uses a simple R script to visualize resource usage of a resource-hungry task."
    }
    parameter_meta {
        resource_log: "Resource usage log file."
        output_pdf_name: "name of the output plot, must end in .pdf"
        plot_title: "Title of the plot"
    }
    input {
        File resource_log
        String output_pdf_name
        String plot_title
    }
    output {
        File plot_pdf = "~{output_pdf_name}"
    }

    command <<<
        set -eux

        /opt/plot.resources.R "~{resource_log}" "~{output_pdf_name}" "~{plot_title}"
    >>>
    runtime {
        disks: "local-disk 10 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-resource-visual:0.1.1"
    }
}
