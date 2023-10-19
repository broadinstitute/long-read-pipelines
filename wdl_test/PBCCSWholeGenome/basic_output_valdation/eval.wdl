version 1.0

workflow eval_workflow {
    input {

        #Floats

        Float aligned_read_length_N50
        Float aligned_num_reads
        Float aligned_frac_bases
        Float aligned_num_bases
        Float aligned_read_length_stdev
        Float average_identity
        Float aligned_est_fold_cov
        Float aligned_read_length_mean
        Float median_identity
        Float aligned_read_length_median

        #Files

        String pbsv_tbi
        String sniffles_vcf
        String clair_gtbi
        String dvp_tbi
        String dvp_g_tbi
        String dvp_vcf
        String clair_vcf
        String pbsv_vcf
        String aligned_pbi
        String aligned_bai
        String dvp_phased_vcf
        #String bed_cov_summary # this ends up being 'null' so not including in array of files for now
        String dvp_phased_tbi
        String clair_tbi
        String clair_gvcf
        String aligned_bam
        String sniffles_tbi
        String dvp_g_vcf

    }

    Array[Float] workflow_out_floats = [
        aligned_read_length_N50,
        aligned_num_reads,
        aligned_frac_bases,
        aligned_num_bases,
        aligned_read_length_stdev,
        average_identity,
        aligned_est_fold_cov,
        aligned_read_length_mean,
        median_identity,
        aligned_read_length_median
    ]
    Array[String] workflow_out_files = [
        pbsv_tbi,
        sniffles_vcf,
        clair_gtbi,
        dvp_tbi,
        dvp_g_tbi,
        dvp_vcf,
        clair_vcf,
        pbsv_vcf,
        aligned_pbi,
        aligned_bai,
        dvp_phased_vcf,
        dvp_phased_tbi,
        clair_tbi,
        clair_gvcf,
        aligned_bam,
        sniffles_tbi,
        dvp_g_vcf
    ]

    String ubuntu_image_tag = "latest"
    String gcloud_slim_image_tag = "slim"

################
    ## Compairing test and expected Floats
    ## Json inputs like "eval_workflow.workflow_out_floats": [{"Left":1.0,"Right":1.0},{"Left":2.2,"Right":3.2}]
    ## With Array[Pair[Float,Float]] workflow_out_floats

#    Array[Boolean] scattered_float_match = []
#    scatter (pair in workflow_out_floats){
#        if (pair.left != pair.right) {
#            Boolean scattered_float_match = false
#            call CheckerWorkflowError{
#                input:
#                    message = "Expected Float "+pair.right+" but got "+pair.left,
#                    image_to_use = ubuntu_image
#            }
#        }
#    }
################

    ## Confirm float does not equal to zero test and expected Floats
    Array[Boolean] scattered_float_match = []
    scatter (in_float in workflow_out_floats){
        if (in_float == 0.0) {
            Boolean scattered_float_match = false
            call CheckerWorkflowError{
                input:
                    message = "Error: Expected a non-zero float but got "+in_float+" .",
                    image_tag = ubuntu_image_tag
            }
        }
    }


    call CheckFileUpdatedDateGCP {
        input:
        file_paths = workflow_out_files,
        image_tag = gcloud_slim_image_tag
    }
}

task CheckFileUpdatedDateGCP {

    meta {
        description: "Checks if file was updated within a specified time (default: 1 day)"
        note: "Specific to GCP"
    }


    input {
        Array[String] file_paths
        Int days_back = 1
        String image_tag
    }

    String image_to_use = "gcr.io/google.com/cloudsdktool/cloud-sdk:" + image_tag

    #FILE_DATE description: get file info | grep the 'Update Time' row | parse the date info | reformat date info

    command <<<
        set -eu pipefail

        GS_BUCKET_PATHS=("~{sep='" "' file_paths}")
        EMPTY_MD5="d41d8cd98f00b204e9800998ecf8427e"
        return_code=0
        echo -e "FilePath\tFileMD5\tValid"

        for GS_FILE in ${GS_BUCKET_PATHS[@]};
        do
          FILE_MD5=`gsutil hash -hm  ${GS_FILE} | grep "md5" | awk -F '[:][\t]+' '{print $2}'`

          if [[ $FILE_MD5 != $EMPTY_MD5 ]] ; then
            echo -e "$GS_FILE\t$FILE_MD5\tTrue"
          else
            echo "ERROR: MD5 for $GS_FILE equals md5sum of an empty file: $EMPTY_MD5" >&2

            echo -e "$GS_FILE\t$FILE_MD5\tFalse"
            return_code=1
          fi
        done

        if [ $return_code == 1 ]; then
            exit 1
        fi

    >>>
    runtime {
        docker: image_to_use
    }
    output {
        File file_date_result = stdout()
    }
}

task CheckerWorkflowError {

    input {
        String message
        String image_tag
    }

    String image_to_use = "marketplace.gcr.io/google/ubuntu2004:" + image_tag

    command <<<
        set -eu pipefail

        echo ~{message}
        exit 1

    >>>
    runtime {
        docker: image_to_use
    }
    output {
        Boolean errmessage = stdout()
    }
}

task ValidFloatOutput {

    input {
        Array[Pair[Float,Float]] workflow_out_floats
        String image_to_use
    }
    command <<<

        echo ~{workflow_out_floats}

    >>>
    runtime {
        docker: image_to_use
    }
    output {
        Boolean comparison_result = read_boolean(stdout())
    }
}

task ValidMd5SumOutput {
    input {
        File data_file
        String expectedMd5sum
    }
    command <<<

        md5sum helloworld.txt | sed "s/|/ /" | awk "{print $1, $8}" | read filemd5

        if [$filemd5 == $expectedMd5sum]
        then
        echo "true"
        else
        echo "false"
        fi
    >>>
    runtime {
        docker: "quay.io/agduncan94/my-md5sum"
    }
    output {
        File comparison_result = stdout()
    }
}
