version 1.0

workflow TerraEnvironmentTesting {
    input {
        String? test_string
    }

    call CheckEnvironmentVariables{input:}

    output {
        File environment_variables = CheckEnvironmentVariables.environment_variables
    }
}

task CheckEnvironmentVariables {
    input {
    }

    command {
        env
        env > environment_variables.txt
    }

    output {
        File environment_variables = "environment_variables.txt"
    }

    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2204:b25bdce835cb"
    }
}