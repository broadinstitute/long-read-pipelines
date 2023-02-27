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

task StreamingBamErrored {
    input {
        File stderr_log
    }
    command <<<
        set -eux
        if [[ -s ~{stderr_log} ]]; then
            echo "Streaming a BAM triggered warnings or errors." \
            && echo "true" > result.txt \
            && cat ~{stderr_log} \
            && exit 1
        fi
    >>>
    output {
        Boolean yes = read_boolean("result.txt")
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

task StreamCountAlignmentRecords {
    meta {
        desciption: "Count the number of alignment records with a particular SAM flag. The BAM is streamed in, so the task may fail potentially."
    }
    input {
        File aligned_bam
        File aligned_bai

        Int? decimal_flag
    }

    parameter_meta {
        decimal_flag: "The SAM flag, in decimal (as opposed to hexadecimal.)"

        aligned_bam: {
            localization_optional: true
        }
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -c \
            ~{true='-f ' false=' ' defined(decimal_flag)}~{select_first([decimal_flag, " "])} \
            ~{aligned_bam} \
            > "count.txt" \
            2>"error.log"
        touch "error.log"
    >>>

    output {
        File stderr_log =  "error.log"
        Int count = read_int("count.txt")
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

task CountAlignmentRecordsByFlag {
    meta {
        desciption: "Count the number of alignment records with a particular SAM flag."
    }
    input {
        File aligned_bam
        File aligned_bai

        Map[String, Int] names_and_decimal_flags

        Int num_local_ssds
    }

    parameter_meta {
        names_and_decimal_flags: "The SAM flags, in decimal (as opposed to hexadecimal.), with their (appropriate) names."
    }

    # Array[Pair[String, Int]] dummy = names_and_decimal_flags
    Int n = length(names_and_decimal_flags)

    Int disk_gb = 375 * num_local_ssds

    command <<<
        set -eux

        # iterate through each requested SAM flag
        two_col_tsv=~{write_map(names_and_decimal_flags)}
        cat "${two_col_tsv}"
        x=$(wc -l "${two_col_tsv}" | awk '{print $1}')
        if [[ ~{n} -ne "${x}" ]]; then
            sed -i -e '$a\' "${two_col_tsv}" # because write_map doesn't have new line at end of file, so we add it explicitly
        fi
        wc -l "${two_col_tsv}"

        # filter and count'
        while IFS=$'\t' read -r -a line
        do
            name="${line[0]}"
            flag="${line[1]}"
            samtools view -c \
                -f "${flag}" \
                ~{aligned_bam} \
                > "asdfxyz_${name}.txt" &
        done < "${two_col_tsv}"
        # primary
        samtools view -c \
            -F 2308 \
            ~{aligned_bam} \
            > "primary_count.txt" &
        # overall
        samtools view -c ~{aligned_bam} > "total_cnt.txt" &
        wait
        total_count=$(cat total_cnt.txt)
        primary_count=$(cat primary_count.txt)

        # format results
        touch result_raw.txt result_pct.txt # have to do this because iterating a Map isn't possible on WDL 1.0
        echo -e "Primary\t${primary_count}" > result_raw.txt
        pct=$(echo "100*${primary_count}/${total_count}" | bc | awk '{ printf("%.0f\n",$1) }')
        echo -e "Primary\t${pct}" > result_pct.txt
        for ff in asdfxyz*txt;
        do
            name=$(echo "${ff}" | sed 's#asdfxyz_##' | sed 's#.txt$##')
            count=$(cat "${ff}")
            echo -e "${name}\t${count}" >> result_raw.txt
            pct=$(echo "100*${count}/${total_count}" | bc -l | awk '{ printf("%.0f\n",$1) }')
            echo -e "${name}\t${pct}" >> result_pct.txt
        done

        awk -F '\t' '{print $2}' result_pct.txt | paste -sd+ | bc > debug_total_pct.txt
    >>>

    output {
        Int total_cnt = read_int("total_cnt.txt")
        Map[String, Int] flag_cnts = read_map("result_pct.txt")
        Map[String, Float] flag_pcts = read_map("result_pct.txt")
        Int summed_percentages = read_int("debug_total_pct.txt")
    }

    runtime {
        cpu:            8
        memory:         "32 GiB"
        disks:          "local-disk " + disk_gb + " LOCAL"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.2"
    }
}
