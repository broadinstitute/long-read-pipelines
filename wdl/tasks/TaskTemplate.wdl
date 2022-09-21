## #NameofTaskWorkflow
##
## Task Workflow Description, e.g.:
## This WDL defines tasks used for ...

version 1.0

task TaskNameOne {
    meta {
        description: "Task description"
    }
    parameter_meta {
        input_one:   "Input description"
        input_two:   "Input description"
    }

    input {
        File input_one
        String input_two
    }

    command {
        set -eux
        echo "Hello World"
    }

    output {
        String output_one = stdout()
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }

}

task TaskNameTwo {
    meta {
        description: "Task description"
    }
    parameter_meta {
        input_one:   "Input description"
        input_two:   "Input description"
    }

    input {
        File input_one
        String input_two
    }

    command {
        set -eux
        echo "Hello World"
    }

    output {
        String output_one = stdout()
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }

}

task TaskNameThree {
    meta {
        description: "Task description"
    }
    parameter_meta {
        input_one:   "Input description"
        input_two:   "Input description"
    }

    input {
        File input_one
        String input_two
    }

    command {
        set -eux
        echo "Hello World"
    }

    output {
        String output_one = stdout()
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }

}
