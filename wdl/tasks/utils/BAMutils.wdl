version 1.0

task GetReadGroupInfo {
    meta {
        desciption:
        "Get some read group information Given a single-readgroup BAM. Will fail if the information isn't present."
    }

    input {
        String uBAM  # not using file as call-caching brings not much benefit

        Array[String] keys
    }

    parameter_meta {
        keys: "A list of requested fields in the RG line, e.g. ID, SM, LB."
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

task CountReadGroupsAndSamples {
    input {
        String  bam
        String? bai
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} | grep "^@RG" > rg_header.txt

        wc -l rg_header.txt | awk '{print $1}' > rg_cnt.txt

        cat rg_header.txt | tr '\t' '\n' | grep -cF "SM:" > sm_cnt.txt
    >>>

    output {
        Int sm_cnt = read_int("sm_cnt.txt")
        Int rg_cnt = read_int("rg_cnt.txt")
    }

    runtime {
        cpu:            2
        memory:         "8 GiB"
        disks:          "local-disk 100 HDD"
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task ReSampleBam {
    meta {
        description: "For resetting the sample name in a BAM file to the requested value"
    }
    input {
        File  bam
        File? bai

        String new_sample_name
        String out_prefix
    }
    parameter_meta {
        bam: "BAM to operate on, aligned or not, but is assumed to contain @RG lines in its header holding SM information, and the BAM is for a single sample."
        bai: "Accompanying bai index."
        new_sample_name: "New sample name desired (assumed to be a valid value)."
        out_prefix: "Prefix for output BAM."
    }

    Boolean index = defined(bai)

    command <<<
        set -eux

        samtools view -H ~{bam} > header.txt
        echo "================"
        grep "^@RG" header.txt
        echo "================"

        grep "@RG" header.txt | tr '\t' '\n' | grep "^SM:" | awk -F ':' '{print $2}' \
        > all_current_sample_names.txt

        while IFS= read -r old_sample_name
        do
            echo "${old_sample_name}"
            sed -i.bak \
                "s/SM:${old_sample_name}/SM:~{new_sample_name}/g" \
                header.txt
        done < all_current_sample_names.txt
        echo "================"
        grep "^@RG" header.txt
        echo "================"

        samtools reheader header.txt ~{bam} > ~{out_prefix}.bam
        if ~{index} ; then samtools index -@1 ~{out_prefix}.bam; fi
    >>>

    output {
        File  reheadered_bam = "~{out_prefix}.bam"
        File? reheadered_bai = "~{out_prefix}.bam.bai"
    }

    runtime {
        cpu:            2
        memory:         "8 GiB"
        disks:          "local-disk 375 LOCAL"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}
