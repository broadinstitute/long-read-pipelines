## #

## #
version 1.0

workflow NameofWorkflow2 {

    meta {
        description:   "Workflow description"
        outputs:   "- `output_one` Output description"
    }

    parameter_meta {
        input_one:   "Input description"
        input_two:   "Input description"
        input_three: "Input description"
    }

    input {
        File input_one
        String input_two

        Int? input_three = 10

    }

    call TaskName {
        input:
            input_one = input_one,
            input_two = input_two
    }

    output {
        Float output_one = TaskName.output_one
    }

}

task TaskName {
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
