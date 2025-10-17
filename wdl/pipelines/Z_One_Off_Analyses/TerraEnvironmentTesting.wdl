version 1.0

workflow TerraEnvironmentTesting {
    input {
    }

    call CheckEnvironmentVariables{input:}

    output {
        File environment_variables = CheckEnvironmentVariables.environment_variables
    }
}

task CheckEnvironmentVariables {


    meta {
        volatile: true
    }

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
        docker: "gcr.io/cloud-marketplace/google/ubuntu2204:b25bdce835cbcde2cf02c94dd37cca3c6b43c6967f5a8a8a8e6a77e4956551a0"
    }
}