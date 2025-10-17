version 1.0

import "../../tasks/Utility/Utils.wdl" as UTILS

workflow TerraEnvironmentTesting {
    input {
    }

    call UTILS.check_terra_env as check_terra_env {}
    
    output {
        Boolean is_running_on_terra       = check_terra_env.is_running_on_terra

        Boolean is_backed_by_gcp          = check_terra_env.is_backed_by_gcp
        Boolean is_running_via_gcp_batch  = check_terra_env.is_running_via_gcp_batch
        Boolean is_running_via_gcp_papiv2 = check_terra_env.is_running_via_gcp_papiv2

        String google_project_id          = check_terra_env.google_project_id

        String user_email                 = check_terra_env.user_email

        String workspace_uuid             = check_terra_env.workspace_uuid
        String workspace_name             = check_terra_env.workspace_name
        String workspace_namespace        = check_terra_env.workspace_namespace
        String workspace_bucket_path      = check_terra_env.workspace_bucket_path

        String method_version             = check_terra_env.method_version
        String method_source              = check_terra_env.method_source
        String method_path                = check_terra_env.method_path

        #String gcp_created_by_metadata = read_string("gcp_created_by_attributes.txt")
        File   gcp_instance_metadata      = check_terra_env.gcp_instance_metadata

        String input_table_name           = check_terra_env.input_table_name
        String input_row_id               = check_terra_env.input_row_id

        String top_level_submission_id    = check_terra_env.top_level_submission_id

        File env_info                     = check_terra_env.env_info
        File gcloud_config_info           = check_terra_env.gcloud_config_info
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
        docker: "gcr.io/cloud-marketplace/google/ubuntu2204@sha256:b25bdce835cbcde2cf02c94dd37cca3c6b43c6967f5a8a8a8e6a77e4956551a0"
    }
}