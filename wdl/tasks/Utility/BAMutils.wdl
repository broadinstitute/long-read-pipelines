version 1.0

task GetReadGroupInfo {
    meta {
        desciption:
        "Get some read group information Given a single-readgroup BAM. Will fail if the information isn't present."
    }

    parameter_meta {
        uBAM: "The input BAM file."
        keys: "A list of requested fields in the RG line, e.g. ID, SM, LB."
    }

    input {
        String uBAM  # not using file as call-caching brings not much benefit

        Array[String] keys
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{uBAM} | grep "^@RG" | tr '\t' '\n' > rh_header.txt

        for attribute in ~{sep=' ' keys}; do
            value=$(grep "^${attribute}" rh_header.txt | awk -F ':' '{print $2}')
            echo -e "${attribute}\t${value}" >> "result.txt"
        done
    >>>

    output {
        Map[String, String] read_group_info = read_map("result.txt")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}
