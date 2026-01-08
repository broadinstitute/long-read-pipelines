version 1.0

import "../../structs/Structs.wdl"

#################################################
# header-only READ operations
#################################################

task GetReadGroupInfo {
    meta {
        desciption:
        "Get some read group information given a single-readgroup BAM. If the requested keys are absent, a null value is assigned in the returned entry."
        warn:
        "If the BAM contains multiple read groups, task will fail."
    }
    parameter_meta {
        bam: { desciption: "The input BAM file.", localization_optional: true }
        keys: "A list of requested fields in the RG line, e.g. ID, SM, LB."
        null_value_representation: "For keys requested that aren't available in the bam's header, this value will be returned."
    }

    input {
        File bam

        Array[String] keys
        String null_value_representation = "None"
    }

    output {
        Map[String, String] read_group_info = read_map("result.tsv")
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} | grep "^@RG" > one_rg_per_line.txt
        num_rgs=$(wc -l one_rg_per_line.txt | awk '{print $1}')
        if [[ ${num_rgs} -gt 1 ]]; then exit 1; fi

        cat one_rg_per_line.txt | tr '\t' '\n' > rh_header.txt

        for attribute in ~{sep=' ' keys}; do
            if grep -q "^${attribute}" rh_header.txt; then
                value=$(grep "^${attribute}" rh_header.txt | awk -F ':' '{print $2}')
            else
                value="~{null_value_representation}"
            fi
            echo -e "${attribute}\t${value}" >> "result.tsv"
        done
    >>>

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 10 HDD"
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}

task GetReadGroupLines {
    meta {
        desciption: "Get the @RG lines in a BAM's header. Will error if there's no read group defined in the header."
    }
    parameter_meta {
        bam: {localization_optional: true}
    }

    input {
        File bam
    }

    output {
        Array[String] read_group_ids = read_lines("rgids.txt")
        Array[String] read_group_lines = read_lines("read_groups.txt")
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        samtools view -H ~{bam} | grep "^@RG" > read_groups.txt

        rm -f rgids.txt
        while IFS= read -r line
        do
            echo "${line}" | tr '\t' '\n' \
                | grep "^ID:" | cut -d ':' -f 2- \
            >> rgids.txt
        done < read_groups.txt
    >>>

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 10 HDD"
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}

task GatherBamMetadata {
    meta {
        description: "Check several metadata of an input BAM (aliged? sort order? etc)"
    }
    parameter_meta {
        bam: { localization_optional: true }
    }

    input {
        File bam
    }

    output {
        Boolean is_aligned = read_boolean("is_mapped.txt")

        Boolean is_sorted = read_boolean("is_sorted.txt")
        String sort_order = read_string("sort_order.txt")
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} > header.txt

        grep -F "@HD" header.txt | tr '\t' '\n' > hd.line.txt
        if grep -q "SO" hd.line.txt;
        then
            echo "true" > "is_sorted.txt"
            grep "SO" hd.line.txt | cut -d ':' -f 2- \
            > "sort_order.txt"
        else
            echo "false" > "is_sorted.txt"
            echo "NA" > "sort_order.txt"
        fi

        # we use two conditions: @SQ lines in header, and at least some mapped reads
        mapped_bool=''
        if grep -q "@SQ" "header.txt";
        then
            export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
            if [[ 1 -le $(samtools view -F 4 ~{bam} | head | wc -l | awk '{print $1}') ]];
            then
                mapped_bool='true'
            else # @SQ lines defined but seemingly no mapped reads
                mapped_bool='unsure' # this will trigger error in WDL later, but it's intentional because we cannot be sure
            fi
        else
            mapped_bool='false'
        fi
        echo "${mapped_bool}" > "is_mapped.txt"
    >>>

    runtime {
        cpu:    1
        memory: "4 GiB"
        disks:  "local-disk 10 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}

task CountReadGroups {
    meta {
        desciption: "Count the number of RG lines in the header of the BAM file."
    }
    parameter_meta {
        bam: { localization_optional: true }
    }

    input {
        File bam
    }

    output {
        Int num_rg = read_int("rg_cnt.txt")
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} | grep -c "^@RG" > "rg_cnt.txt"
    >>>

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 10 HDD"
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}

task InferSampleName {
    meta {
        description: "Infer sample name encoded on the @RG line of the header section."
        warn: "Fails if multiple values found, or if SM is the specified illegal value."
    }
    parameter_meta {
        bam: { localization_optional: true }
    }

    input {
        File bam
        File? bai
        String illegal_value = "unnamedsample"
    }

    output {
        String sample_name = read_string("sample.names.txt")
    }

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} > header.txt
        if ! grep -q '^@RG' header.txt; then echo "No read group line found!" && exit 1; fi

        grep '^@RG' header.txt | tr '\t' '\n' | grep '^SM:' | sed 's/SM://g' | sort | uniq > sample.names.txt
        if [[ $(wc -l sample.names.txt) -gt 1 ]]; then echo "Multiple sample names found!" && exit 1; fi
        if grep -iq "~{illegal_value}" sample.names.txt; then echo "Sample name found to be illegal!" && exit 1; fi
    >>>

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        preemptible:    2
        maxRetries:     1
        docker:         "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}

#################################################
# statistics
#################################################

task ValidateSamFile {
    meta {
        desciption: "Call GATK/Picard ValidateSamFile to validate input BAM: https://bit.ly/3JMutxp."
    }
    parameter_meta {
        validation_mode: "Desired valiation mode; see Picard documentation for the supproted values."
        disk_type: "Type of disk to use for the computation; SSD for persistent SSD disks, LOCAL for local SSDs."
        bam: {
            localization_optional : true
        }
    }

    input {
        File bam
        String validation_mode = "SUMMARY"

        Array[String] validation_errs_to_ignore = ["INVALID_TAG_NM",  # for the purpose we currently have, NM and CIGAR don't matter, and longreads have no mates
                                                    "MISSING_TAG_NM",
                                                    "INVALID_CIGAR",
                                                    "ADJACENT_INDEL_IN_CIGAR",
                                                    "CIGAR_MAPS_OFF_REFERENCE",
                                                    "MISMATCH_MATE_CIGAR_STRING",
                                                    "MATE_CIGAR_STRING_INVALID_PRESENCE",
                                                    "MATE_NOT_FOUND",
                                                    "INVALID_MAPPING_QUALITY",
                                                    "INVALID_FLAG_MATE_UNMAPPED",
                                                    "MISMATCH_FLAG_MATE_UNMAPPED",
                                                    "INVALID_FLAG_MATE_NEG_STRAND",
                                                    "MISMATCH_FLAG_MATE_NEG_STRAND",
                                                    "INVALID_MATE_REF_INDEX",
                                                    "MISMATCH_MATE_REF_INDEX",
                                                    "MISMATCH_MATE_ALIGNMENT_START",
                                                    "MATE_FIELD_MISMATCH",
                                                    "PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND"
                                                   ]

        String disk_type = "SSD"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(bam, "GiB")) + 50
    String output_basename = basename(basename(bam, ".bam"), ".cram")
    String output_name = "${output_basename}_${validation_mode}.txt"

    String base = basename(bam, ".bam")
    String local_bam = "/mnt/disks/cromwell_root/~{base}.bam"

    command <<<
        set -eux

        time gcloud storage cp ~{bam} ~{local_bam}

        gatk ValidateSamFile \
            --INPUT ~{local_bam} \
            --OUTPUT ~{output_name} \
            --MODE ~{validation_mode} \
            ~{true="--IGNORE " false="" 0<length(validation_errs_to_ignore)} \
            ~{sep=" --IGNORE " validation_errs_to_ignore}
    >>>

    output {
        File validation_report = "${output_name}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-custom-gatk:4.4.0.0-samtools1.18"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SamtoolsFlagStats {
    meta {
        description: "Collect SAM flag stats of an aligned BAM"
    }
    parameter_meta {
        output_format: "argument passed on to '-O' of `samtools flagstats` "
    }

    input {
        File bam
        String output_format = "tsv"
        String disk_type = "SSD"
        RuntimeAttr? runtime_attr_override
    }

    output {
        File flag_stats = "~{output_name}"
    }

    String base = basename(bam, ".bam")

    Map[String, String] reformat_user_input = {'JSON': 'json',       'json': 'json',
                                               "TSV": 'tsv',         'tsv': 'tsv',
                                               'DEFAULT': 'defalt',  'default': 'defalt'}
    String o_f = reformat_user_input[output_format]
    Map[String, String] reformat_output_format = {'default': 'txt', 'json': 'json', 'tsv': 'tsv'}
    String o_ext = reformat_output_format[o_f]

    String output_name = "~{base}.flag_stats.~{o_ext}"

    command <<<
    set -euxo pipefail

        time \
        samtools flagstat \
            -O ~{o_f} \
            ~{bam} \
        > "~{output_name}"
    >>>

    #########################
    Int disk_size = 10 + ceil(size(bam, "GiB"))

    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ParseFlagStatsJson {
    meta {
        description: "Parse the output from `samtools flatstats -O json`"
    }

    parameter_meta {
        sam_flag_stats_json: "JSON output from samtools flatstats"
        replace_none_with: "value to use to replace 'None' in the json value fields"
    }
    input {
        File sam_flag_stats_json
        Float replace_none_with = 0
    }
    output {
        Map[String, Float] qc_pass_reads_SAM_flag_stats = read_map("qcPass.stats.tsv")
        Map[String, Float] qc_fail_reads_SAM_flag_stats = read_map("qcFail.stats.tsv")
    }

    command <<<
        set -euxo pipefail

        python <<CODE
        import json
        with open("~{sam_flag_stats_json}") as inf:
            dd = json.load(inf)
        with open('qcPass.stats.tsv', 'w') as outf:
            for k, v in dd['QC-passed reads'].items():
                vv = "~{replace_none_with}" if v is None else v
                outf.write(f'{k}\t{vv}\n')
        with open('qcFail.stats.tsv', 'w') as outf:
            for k, v in dd['QC-failed reads'].items():
                vv = "~{replace_none_with}" if v is None else v
                outf.write(f'{k}\t{vv}\n')
        CODE
    >>>

    runtime {
        disks: "local-disk 10 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/python:3.9.18-slim-bullseye"
    }
}

task CountMethylCallReads {
    meta {
        desciption: "Count the numbers of records in the bam with and without the ML & MM tags"
    }

    parameter_meta {
        disk_type: "must be one of [HDD, SSD, LOCAL]. SSD is recommended"

        raw_count: "number of records in input BAM"
        bean_count: "number of records in input BAM that has both MM & ML tags"

        non_2304_count: "number of records in input BAM that are neither 256 (secondary) nor 2048 (supplementary)"
        non_2304_bean_count: "number of records in input BAM that are neither 256 (secondary) nor 2048 (supplementary), and have both MM & ML tags"
    }
    input {
        File  bam
        File? bai
        String disk_type
    }
    output {
        Int raw_count  = read_int("raw_count.txt")
        Int bean_count = read_int("bean_count.txt")
        Int non_2304_count = read_int("non_2304_count.txt")
        Int non_2304_bean_count = read_int("non_2304_bean_count.txt")
    }

    command <<<
    set -euxo pipefail

        samtools view -@1 -c         ~{bam} >      raw_count.txt &
        samtools view -@1 -c -F 2304 ~{bam} > non_2304_count.txt &

        samtools view -h         --tag "MM" ~{bam} | samtools view -c --tag "ML" >          bean_count.txt &
        samtools view -h -F 2304 --tag "MM" ~{bam} | samtools view -c --tag "ML" > non_2304_bean_count.txt &

        wait

        tail ./*_count.txt
    >>>
    Int bam_sz = ceil(size(bam, "GiB"))
    Int local_ssd_sz = if bam_sz > 300 then 750 else 375
    Int pd_sz = 50 + bam_sz
    Int disk_size = if "LOCAL" == disk_type then local_ssd_sz else pd_sz
    runtime {
        cpu:            10
        memory:         "20 GiB"
        disks:          "local-disk ~{disk_size} ~{disk_type}"
        preemptible:    1
        maxRetries:     0
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}

task CountAlignmentRecords {
    meta {
        desciption:
        "Count the number of alignment records with a particular SAM flag."
    }

    parameter_meta {
        aligned_bam: {
            localization_optional: true
        }
        localize_bam: "If false, the BAM is streamed in from the bucket instead of localized, but the operation is subject to network instabilities and may also timeout"

        decimal_flag: "Only include records having the SAM flag; in decimal (as opposed to hexadecimal); when not provided, all records are included."
        inverse_flag: "If true, filter away records that has the requested decimal_flag; no effect when decimal_flag is not provided"
    }

    input {
        File aligned_bam
        File aligned_bai

        Boolean localize_bam = false

        Int? decimal_flag
        Boolean inverse_flag = false
    }

    output {
        File stderr_log =  "error.log"
        Int count = read_int("count.txt")
    }

    command <<<
        set -eux

        touch "error.log"

        if ~{defined(decimal_flag)}; then
            if ~{inverse_flag}; then
                filter_op='-F '
            else
                filter_op='-f '
            fi
        else
            filter_op=' '
        fi
        if ~{localize_bam}; then
            time \
            gcloud storage cp ~{aligned_bam} input_bam.bam
            mv ~{aligned_bai} input_bam.bam.bai

            samtools view -c \
                "${filter_op}" ~{select_first([decimal_flag, " "])} \
                input_bam.bam \
                > "count.txt"
        else
            export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
            samtools view -c \
                "${filter_op}" ~{select_first([decimal_flag, " "])} \
                ~{aligned_bam} \
                > "count.txt" \
                2>"error.log"
        fi
    >>>

    Int disk_size = 20 + ceil(size(aligned_bam, "GiB"))
    String disk_type = if localize_bam then "SSD" else "HDD"
    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk ~{disk_size} ~{disk_type}"
        preemptible:    1
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}

task StreamingBamErrored {
    meta {
        desciption:
        "A helper task that reads the error log from task to report if streaming BAM from GCS bucket failed."
    }
    parameter_meta {
        stderr_log: "The stderr log output from task StreamCountAlignmentRecords"
        yes: "if true, then the streaming read operation failed and the output of that task is likely corrupted."
    }
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
        else
            echo "false" > result.txt
        fi
    >>>
    output {
        Boolean yes = read_boolean("result.txt")
    }

    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task CountAlignmentRecordsByFlag {
    meta {
        desciption:
        "Count the number of alignment records with particular SAM flags; the difference to CountAlignmentRecords is this allows specifying multiple SAM flags"
    }

    parameter_meta {
        names_and_decimal_flags: "The SAM flags, in decimal (as opposed to hexadecimal.), with their (appropriate) names."

        aligned_bam: {
            localization_optional: true
        }
    }

    input {
        File aligned_bam
        File aligned_bai

        Map[String, Int] names_and_decimal_flags

        Int num_local_ssds
    }

    # Int n = length(names_and_decimal_flags)

    String base = basename(aligned_bam, ".bam")
    String local_bam = "/mnt/disks/cromwell_root/~{base}.bam"
    command <<<
        set -eux

        two_col_tsv=~{write_map(names_and_decimal_flags)}
        cat "${two_col_tsv}"
        # x=$(wc -l "${two_col_tsv}" | awk '{print $1}')
        # if [[ 3 -ne "${x}" ]]; then ## 3 here is used in place of n to avoid validation error (yes, validation false positive)
        #     sed -i -e '$a\' "${two_col_tsv}"
        # fi
        wc -l "${two_col_tsv}"
        # because some Cromwell versions' stdlib function write_map() doesn't have new line at end of file, so we add it explicitly
        if [[ $(tail -c1 "${two_col_tsv}" | wc -l) -eq 0 ]]; then
            sed -i -e '$a\' "${two_col_tsv}"
        fi
        # '
        wc -l "${two_col_tsv}"

        time \
            gcloud storage cp ~{aligned_bam} ~{local_bam}
        mv ~{aligned_bai} ~{local_bam}.bai

        # iterate through each requested SAM flag
        while IFS=$'\t' read -r -a line
        do
            name="${line[0]}"
            flag="${line[1]}"
            samtools view -c \
                -f "${flag}" \
                ~{local_bam} \
                > "asdfxyz_${name}.txt" &
        done < "${two_col_tsv}"

        # overall
        samtools view -c ~{local_bam} \
        > "total_cnt.txt" &
        # primary (2308=4+256+2048)
        samtools view -c \
            -F 2308 \
            ~{local_bam} \
        > "primary_count.txt" &
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

    Int local_ssd_sz = if size(aligned_bam, "GiB") > 300 then 750 else 375
    runtime {
        cpu:            8
        memory:         "32 GiB"
        disks:          "local-disk " + local_ssd_sz + " LOCAL"
        preemptible:    1
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}

task GetDuplicateReadnamesInQnameSortedBam {
    meta {
        desciption: "Get read names from a queryname-sorted bam, where such reads are duplicate records"
    }
    parameter_meta {
        qns_bam: {
            localization_optional: true
        }
    }
    input {
        File qns_bam
    }

    output {
        File dup_names_txt = "dup_read_names.txt"
        Boolean result_may_be_corrupted = read_boolean("samtools.failed.txt")
    }

    command <<<
        # the way this works is the following:
        # 0) relying on the re-auth.sh script to export the credentials
        # 1) perform the remote sam-view subsetting in the background
        # 2) listen to the PID of the background process, while re-auth every 1200 seconds
        source /opt/re-auth.sh
        set -euxo pipefail

        # assumption
        sort_order=$(samtools view -H ~{qns_bam} | grep "^@HD" | tr '\t' '\n' | grep "^SO:" | awk -F ':' '{print $2}')
        if [[ "queryname" != "${sort_order}"  ]]; then echo -e "Sort order ${sort_oder} isn't the expected 'queryname'." && exit 1; fi

        # remote grab read names
        echo "false" > samtools.failed.txt
        samtools view ~{qns_bam} \
        | awk -F '\t' '{print $1}' \
        | uniq -d  \
        > "dup_read_names.txt" \
        || { echo "true" > samtools.failed.txt; exit 77; } &
        pid=$!

        set +e
        count=1
        while true; do
            sleep 1200 && date && source /opt/re-auth.sh
            if [[ ${count} -gt 2 ]]; then exit 0; fi
            if ! pgrep -x -P $pid; then exit 0; fi
            count=$(( count+1 ))
        done
    >>>

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 10 HDD"
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}

#################################################
# light transformations (essentially reheader operations)
#################################################

task ResetSamplename {
    meta {
        desciption:
        "Reset the SM entry in the input bam's readgroup lines."
    }

    parameter_meta {
        bam: {localization_optional: true}
    }

    input {
        File bam
        File? bai
        String sample_name
    }

    output {
        File  reheadered_bam = "~{out_prefix}.bam"
        File? reheadered_bai = "~{out_prefix}.bam.bai"
    }

    String prefix = basename(bam, ".bam")
    String local_bam = "/mnt/disks/cromwell_root/~{prefix}.bam"
    String out_prefix = "~{prefix}.ResetSamplename"
    command <<<
        set -euxo pipefail

        time \
        gcloud storage cp ~{bam} ~{local_bam}
        if ~{defined(bai)}; then touch ~{bai}; mv ~{bai} ~{local_bam}.bai; fi

        ###### cleanup the header
        samtools view --no-PG -H ~{local_bam} > header.txt
        grep -v "^@SQ" header.txt

        # fix SM in the RG lines
        grep "^@RG" header.txt > rg_lines.txt
        if ! grep -qF "SM:" rg_lines.txt; then
            sed -i "s/$/SM:tbd/" rg_lines.txt
        fi
        awk -v sm="~{sample_name}" -F '\t' 'BEGIN {OFS="\t"} { for (i=1; i<=NF; ++i) { if ($i ~ "SM:") $i="SM:"sm } print}' \
            rg_lines.txt \
            > fixed_rg_lines.txt
        cat fixed_rg_lines.txt

        # paste things back
        grep -v "^@RG" header.txt > otherlines.txt
        cat otherlines.txt fixed_rg_lines.txt > fixed_header.txt

        ###### samtools reheader
        time \
        samtools reheader fixed_header.txt ~{local_bam} \
        > "~{out_prefix}.bam"
        if ~{defined(bai)}; then
            time \
            samtools index -@1 "~{out_prefix}.bam"
        fi
    >>>

    Int disk_size = 50 + 2 * ceil(size(bam, "GiB"))
    runtime {
        cpu:            2
        memory:         "8 GiB"
        disks:          "local-disk ~{disk_size} SSD"
        preemptible:    2
        maxRetries:     1
        docker:         "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}

#################################################
# intensive transformations -- filter
#################################################

task FilterBamByLen {
    meta {
        desciption: "Filter a BAM by sequence length, and count the yileld if so requested"
        warn: "It's assumed that for aligned BAM, alignment was done without hard clipping turned on. If this assumption isn't met, the resulting BAM may be corrupt."
    }

    parameter_meta {
        len_threshold_inclusive: "Reads longer than or equal to this length will be included."
        bam : { localization_optional: true }
    }

    input {
        File bam
        File? bai
        Int len_threshold_inclusive

        Boolean compute_yield = false

        String disk_type = "HDD"
        RuntimeAttr? runtime_attr_override
    }

    output {
        File  fBAM = "~{out_prefx}.bam"
        File? fBAI = "~{out_prefx}.bam.bai"

        Float?    total_yield = read_float("all.yield.txt")
        Float? filtered_yield = read_float("filtered.yield.txt")
    }

    String base = basename(bam, ".bam")
    String out_prefx = "~{base}.RL_ge_~{len_threshold_inclusive}"

    Boolean is_aligned = defined(bai)

    String local_bam = "/mnt/disks/cromwell_root/~{base}.bam"

    command <<<
        set -euxo pipefail

        time \
        gcloud storage cp ~{bam} ~{local_bam}
        if ~{defined(bai)}; then mv ~{bai} "~{local_bam}.bai"; touch "~{local_bam}.bai"; fi

        # get total yield in the background
        if ~{compute_yield}; then
            # simply get the length of the sequence, excluding 2304 reads
            samtools view -@1 \
                ~{true='-F2304' false=' ' is_aligned} \
                ~{local_bam} \
            | awk -F '\t' '{print length($10)}' \
            > all.read.lengths.txt &
        fi

        if ~{is_aligned} ; then
            # note here that 2048 and 256 reads are not longer than the primary record,
            # so if the primary is already shorter than the threshold, they should be excluded too
            # on the other hand, it is possible the 2048 records are shorter that the actual query when hardclipping is turned on
            # hence we requre the input doesn't have hardclipping turned on
            samtools view -@1 -h \
                --write-index \
                -e "length(seq)>=~{len_threshold_inclusive}" \
                -o "~{out_prefx}.bam##idx##~{out_prefx}.bam.bai" \
                ~{local_bam}
        else
            samtools view -@1 -h \
                -e "length(seq)>=~{len_threshold_inclusive}" \
                -o "~{out_prefx}.bam" \
                ~{local_bam}
        fi

        if ~{compute_yield}; then
            samtools view -@1 \
                ~{true='-F2304' false=' ' is_aligned} \
                "~{out_prefx}.bam" \
            | awk -F '\t' '{print length($10)}' \
            > filtered.read.lengths.txt

            # see https://stackoverflow.com/questions/450799/shell-command-to-sum-integers-one-per-line
            awk '{s+=$1} END {printf "%.0f", s}' filtered.read.lengths.txt \
            > filtered.yield.txt

            wait # make sure total yield gather in the background is done
            awk '{s+=$1} END {printf "%.0f", s}' all.read.lengths.txt \
            > all.yield.txt
        fi
    >>>
    ###################
    Int disk_size = 20 + 2 * ceil(size(bam, "GiB"))

    RuntimeAttr default_attr = object {
        cpu_cores:             4,
        mem_gb:                16,
        disk_gb:               disk_size,
        preemptible_tries:     1,
        max_retries:           1,
        docker:                "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                   select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:           select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:            select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                select_first([runtime_attr.docker, default_attr.docker])
    }
}

task GatherReadsWithoutMethylCalls {
    meta {
        desciption: "Collect records in the bam without the ML & MM tags"
    }
    parameter_meta {
        disk_type: "must be one of [HDD, SSD, LOCAL]. SSD is recommended"
    }

    input {
        File  bam
        File? bai
        String disk_type
    }

    output {
        File no_ml_reads = "~{p}.no_ML.bam"
        File no_mm_reads = "~{p}.no_MM.bam"

        File names_missing_only_one_tag = "missing_only_one_tag.read_names.txt"
        File names_missing_both_tags    = "no_mm_and_ml.read_names.txt"
    }

    String p = basename(bam, ".bam")

    command <<<
    set -euxo pipefail

        export LC_ALL=C  # attempt to make grep faster
        samtools view -@1 -h ~{bam} \
            | grep -vF "ML:B:C" \
            | samtools view -@1 -bh \
            -o "~{p}.no_ML.bam" &

        samtools view -@1 -h ~{bam} \
            | grep -vF "MM:Z:" \
            | samtools view -@1 -bh \
            -o "~{p}.no_MM.bam" &

        wait

        ##########
        samtools view -@1 "~{p}.no_ML.bam" | awk -F '\t' '{print $1}' | sort > no_ml.txt &
        samtools view -@1 "~{p}.no_MM.bam" | awk -F '\t' '{print $1}' | sort > no_mm.txt &
        wait

        comm -3 \
            no_ml.txt \
            no_mm.txt \
        > "missing_only_one_tag.read_names.txt"
        comm -12 \
            no_ml.txt \
            no_mm.txt \
        > "no_mm_and_ml.read_names.txt"
    >>>

    # here, we are a little "brave", given that the task generates some bams
    # however, if the output bams are large enough, then the input definitely has issues (too many without MM/ML)
    # hence the run will fail and we'll be warned by an OoD (or PAPI 10)
    Int bam_sz = ceil(size(bam, "GiB"))
    Int local_ssd_sz = if bam_sz > 300 then 750 else 375
    Int pd_sz = 50 + bam_sz
    Int disk_size = if "LOCAL" == disk_type then local_ssd_sz else pd_sz
    runtime {
        cpu:            10
        memory:         "20 GiB"
        disks:          "local-disk ~{disk_size} ~{disk_type}"
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}

task SubsetBamToLocusLocal {
    meta {
        description: "For subsetting a BAM stored in GCS, expliciting localizing the BAM"
        note: "This is intended as last resort when streaming from buckets fails"
    }

    parameter_meta {
        interval_list_file:  "a Picard-style interval list file to subset reads with"
        interval_id:         "an ID string for representing the intervals in the interval list file"
        prefix: "prefix for output bam and bai file names"
        bam: {localization_optional: true}
    }

    input {
        File bam
        File bai

        File interval_list_file
        String interval_id
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Array[String] intervals = read_lines(interval_list_file)

    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    String subset_prefix = prefix + "." + interval_id

    String local_bam = "/mnt/disks/cromwell_root/~{basename(bam)}"

    command <<<
        set -euxo pipefail

        time gcloud storage cp ~{bam} ~{local_bam}
        mv ~{bai} "~{local_bam}.bai" && touch "~{local_bam}.bai"

        # see man page for what '-M' means
        samtools view \
            -bhX \
            -M \
            -@ 1 \
            --write-index \
            -o "~{subset_prefix}.bam##idx##~{subset_prefix}.bam.bai" \
            ~{local_bam} "~{local_bam}.bai" \
            ~{sep=" " intervals}
    >>>

    output {
        File subset_bam = "~{subset_prefix}.bam"
        File subset_bai = "~{subset_prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"  # expensive, but much faster
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task DeduplicateQuerynameSortedBam {
    meta {
        desciption: "De-duplicate a queryname sorted bam. The queryname sort can be done either in natural order, or ascii order."
    }
    parameter_meta {
        qnorder_bam: {
            desciption: "queryname sorted BAM",
            localization_optional: true
        }
    }
    input {
        File qnorder_bam
        RuntimeAttr? runtime_attr_override
    }
    output {
        File dedup_bam = "~{base}.dedup.bam"
        File dup_read_names = "duplicated.readnames.txt"
    }

    String base = basename(qnorder_bam, ".bam")
    String local_bam = "/mnt/disks/cromwell_root/~{base}.bam"

    Int disk_size = 3 * ceil(size(qnorder_bam, "GB"))

    command <<<
        set -eux

        time gcloud storage cp ~{qnorder_bam} ~{local_bam}

        # if no duplicate at all, why bother
        time samtools view ~{local_bam} | awk -F '\t' '{print $1}' | sort | uniq -d > duplicated.readnames.txt
        touch duplicated.readnames.txt
        cat duplicated.readnames.txt
        cnt=$(wc -l duplicated.readnames.txt | awk '{print $1}')
        if [[ ${cnt} -eq 0 ]]; then
            echo "No duplicates found in the unmapped reads"
            mv ~{local_bam} "~{base}.dedup.bam"
        else
            time \
            python3 /opt/remove_duplicate_ont_namesorted_unaligned.py \
                -p "~{base}.dedup" \
                -q "duplicated.readnames.bypython.txt" \
                ~{local_bam}

            cat "duplicated.readnames.bypython.txt"

            diff <(sort duplicated.readnames.txt) <(sort duplicated.readnames.bypython.txt)
        fi
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-bam-dedup:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

#################################################
# intensive transformations -- map
#################################################

task BamToFastq {
    meta {
        description : "Convert a long reads BAM file to a fastq file."
        warn: "Please do not include 'RG' in tags_to_preserve, as that's automatically saved"
    }

    parameter_meta {
        bam: {localization_optional: true}
        prefix: "Prefix for the output fastq file."

        save_all_tags:
        "if true, saves all SAM tags to the FASTQ output; cannot set this true while also specifying tags_to_preserve "
        tags_to_preserve:
        "custom list of tags to preserve; please do not include 'RG' in tags_to_preserve, as that's automatically preserved"

        disk_type: "type of disk to use"
    }

    input {
        File bam
        String prefix

        Boolean save_all_tags = false
        Array[String] tags_to_preserve = []

        String disk_type = "SSD"
        RuntimeAttr? runtime_attr_override
    }

    output {
        File reads_fq = "~{prefix}.fq.gz"
    }

    Boolean custom_tags_to_preserve = 0<length(tags_to_preserve)

    String base = basename(bam)
    String local_bam = "/mnt/disks/cromwell_root/~{base}"
    command <<<
        set -euxo pipefail

        # some checks on inputs
        if ~{custom_tags_to_preserve} && ~{save_all_tags} ; then
            echo "cannot ask to save all tags and yet also ask to save a custom list of tags" && exit 1;
        fi
        for tag in ~{sep=' ' tags_to_preserve}; do
            if [[ ${tag} == "RG" ]]; then
                echo "RG tag is automatically saved for you, no need to save it" && exit 1
            fi
        done

        # localize
        time \
        gcloud storage cp ~{bam} ~{local_bam}


        # when saving all tags, the list can be empty as instructed by samtools doc
        time \
        samtools fastq \
            -@1 \
            -t \
            -0 ~{prefix}.fq.gz \
            ~{local_bam}

        # also using pigz to enable parallel compression
        time \
        samtools fastq \
            ~{true='-T  ' false =' ' save_all_tags} \
            ~{true='-T ' false =' ' custom_tags_to_preserve} ~{sep=',' tags_to_preserve} \
            ~{local_bam} \
        | pigz \
        > ~{prefix}.fq.gz
    >>>

    #########################
    Int disk_size = 10 + 3 * ceil(size(bam, "GiB"))

    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task GetPileup {
    meta {
        desciption:
        "Get pileup information with `samtools mpileup`. Current cmdline options are '-a -s -q 1 [-E|-B]' "
        warn:
        "Do not run this on a large BAM, i.e. pre-filter BAM before running this task. Also see if task BamToRelevantPileup is what you need."
    }

    parameter_meta {
        bam: {localization_optional: true}
        disable_baq: "User choice to diable BAQ computation or not (see doc for samtools mpileup for detail)"
    }

    input {
        File bam
        File bai
        Boolean disable_baq
        String prefix
        File ref_fasta
    }

    output {
        File pileup = "~{prefix}.mpileup"
    }

    String baq_option = if disable_baq then '-B' else '-E'

    String base = basename(bam)
    String local_bam = "/mnt/disks/cromwell_root/~{base}"

    command <<<
        set -euxo pipefail

        time gcloud storage cp ~{bam} ~{local_bam}

        samtools mpileup \
            ~{baq_option} \
            -a \
            -s \
            -q 1 \
            -f ~{ref_fasta} \
            -o ~{prefix}.mpileup \
            ~{local_bam}
    >>>

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}

task BamToRelevantPileup {
    meta {
        desciption:
        "Chop up a GRCh38 BAM by chromosome and further subset into requested genotyping sites; then convert to pileup format. See also task GetPileup."
        note:
        "This task may fail due to some strange samtools failures reading from NAS storage (including cloud PDs that aren't local SSDs). So we included inputs and outputs to guard against that."
    }

    parameter_meta {
        bam: {localization_optional: true}
        bed: "sites where pileup is needed"

        max_retries: "Because of the strange samtools failures reading from NAS storage, we should make multiple attempts to get away from the trasient errors. If after the max retries, we still get those failures, this task will fail."

        pileup_stderr: "stderr output from the samtools mpileup commands, they should ALL be 1-liners."
    }
    input {
        File bam
        File bai
        File bed
        File ref_fasta
        Boolean disable_baq

        String disk_type = "SSD"
        Int max_retries = 1
    }
    output {
        File pileups = "pileup.mpileup"
        Array[File] pileup_stderr = glob("*.mpileup.err")
    }

    String baq_option = if disable_baq then '-B' else '-E'

    String base = basename(bam)
    String local_bam = "/mnt/disks/cromwell_root/~{base}"

    command <<<
        set -euxo pipefail

        time \
        gcloud storage cp ~{bam} ~{local_bam}
        mv ~{bai} "~{local_bam}.bai"

        # generate bed for parallel conversion
        set +e
        for i in `seq 1 22`;
        do
            grep -w "chr${i}" ~{bed} > "chr${i}.bed";
        done
        grep -w "chrX" ~{bed} > "chrX.bed"
        grep -w "chrY" ~{bed} > "chrY.bed"
        set -e
        rm ~{bed}

        # parallel conversion
        cnt=0
        for bed in $(ls chr*.bed | sort -V); do

            if [[ ! -s ${bed} ]] ; then rm "${bed}" && continue; fi

            bash /opt/convert.2.pileup.sh \
                ~{local_bam} ~{ref_fasta} ~{baq_option} \
                ${bed} \
                &

            cnt=$((cnt + 1))
            if [[ $cnt -eq ~{cores} ]]; then wait; cnt=0; fi
        done
        wait
        ls -lh

        # here we use a trick, that if any of the stderr file is large, the conversion must have failed
        mpileup_stderr_sz=$(du -c *.mpileup.err | tail -n1 | awk '{print $1}')
        if [[ "${mpileup_stderr_sz}" -gt 1024 ]]; then
            du -c *.mpileup.err
            echo "some chromosome failed to be converted to pileup"
            exit 1
        fi

        rm -f chr*bam chr*bai
        cat *.mpileup > pileup.mpileup
    >>>

    Int cores = 12
    Int memory = 4 + cores
    Int local_ssd_sz = if size(bam, "GiB") > 150 then 750 else 375
    Int pd_sz = 20 + 4 * ceil(size(bam, "GiB"))
    Int min_pd_sz = 300
    Int picked_pd_sz = if pd_sz < min_pd_sz then min_pd_sz else pd_sz
    Int disk_size = if "LOCAL" == disk_type then local_ssd_sz else picked_pd_sz

    runtime {
        cpu:            "~{cores}"
        memory:         "~{memory} GiB"
        disks:          "local-disk ~{disk_size} ~{disk_type}"
        preemptible:    1
        maxRetries:     max_retries
        docker:         "us.gcr.io/broad-dsp-lrma/lr-bam-pileup:0.1.3"
    }
}

task SamtoolsReset {
    meta {
        description: "Use samtools reset to drop alignment information from the input bam."
    }

    parameter_meta {
        bam: {
            desciption: "aligned BAM to operate on",
            localization_optional: true
        }
        addtional_tags_to_drop: "tags in each alignment record to be dropped; usually these are tags produced by the mapper/aligner that generated the original alignment"
    }
    input {
        File bam
        # these are known mm2 tags and pbmm2 tags
        Array[String] addtional_tags_to_drop = ['cg', 'cm', 'cs',
                                                'de', 'dv',
                                                'ms',
                                                'nn',
                                                'rl',
                                                's1', 's2',
                                                'tp', 'ts',
                                                'mc', 'mg', 'mi', 'rm']

        Int? num_ssds
        RuntimeAttr? runtime_attr_override
    }

    output {
        File res = "~{prefix}.unaligned.bam"
        File original_sam_flag_stats = "orignal.SAM-flag.stats.txt"
    }

    Array[String] std_tags_to_drop = ['MD', 'NM', 'AS', 'SA', 'XA']
    Array[String] tags_to_drop = flatten([std_tags_to_drop, addtional_tags_to_drop])

    String prefix = basename(bam, ".bam")

    Int disk_size = if defined(num_ssds) then 375*select_first([num_ssds]) else 1+10*ceil(size([bam], "GB"))
    String disk_type = if defined(num_ssds) then " LOCAL" else " SSD"

    String base = basename(bam, ".bam")
    String local_bam = "/mnt/disks/cromwell_root/~{base}.bam"

    command <<<
        set -eux

        time gcloud storage cp ~{bam} ~{local_bam}

        samtools view -@1 ~{local_bam} | grep -v "^@" | awk -F '\t' '{print $2}' | sort | uniq -c > orignal.SAM-flag.stats.txt &

        samtools reset -@3 \
            --remove-tag ~{sep=',' tags_to_drop} \
            -o ~{prefix}.unaligned.bam \
            ~{local_bam}
        wait
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + disk_type
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task QuerynameSortBamWithSamtools {
    meta {
        description: "queryname-sort a BAM with samtools. WARNING: see https://github.com/samtools/samtools/issues/1500 if you should use samtools"
    }

    parameter_meta {
        bam: "input BAM"
        qnsort_bam: "output BAM sorted by query name"
        num_ssds: "Number of local SSDs to use; if not provided, will use SSD persistent disks"

        multi_record_queries: "File holding names of queries that has multiple records in the output"
    }
    input {
        File bam
        Int? num_ssds
        RuntimeAttr? runtime_attr_override
    }

    output {
        File qnsort_bam = "~{prefix}.qname-sorted.bam"
        File multi_record_queries = "multi_record_queries.txt"
    }

    String prefix = basename(bam, ".bam")

    Int disk_size = if defined(num_ssds) then 375*select_first([num_ssds]) else 1+4*ceil(size([bam], "GB"))
    String disk_type = if defined(num_ssds) then " LOCAL" else " SSD"

    command <<<

        echo "don't use me yet; see if your version of samtools has this ticket resolved https://github.com/samtools/samtools/issues/1500"; exit 1

        set -eux

        samtools view -H ~{bam} | grep "^@HD" | tr '\t' '\n' > hd.line.txt
        if grep -q 'SO:queryname' hd.line.txt;
        then
            echo "already sorted"
            mv ~{bam} "~{prefix}.qname-sorted.bam"
            exit 0
        fi
        samtools sort -@3 -m2G \
            -N \
            -O BAM \
            ~{bam} \
        > "~{prefix}.qname-sorted.bam"

        touch multi_record_queries.txt
        samtools view "~{prefix}.qname-sorted.bam" | awk '{print $1}' | uniq -d > multi_record_queries.txt
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + disk_type
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task QuerynameSortBamWithPicard {
    meta {
        desciption: "See https://github.com/samtools/samtools/issues/1500 why we aren't using samtools. Note that this task is disk-space hungry."
    }

    parameter_meta {
        bam: "input BAM"
        qnsort_bam: "output BAM sorted by query name"
        num_ssds: "Number of local SSDs to use; if not provided, will use SSD persistent disks (instead of local SSDs)"
    }
    input {
        File bam
        Int? num_ssds
        RuntimeAttr? runtime_attr_override
    }

    output {
        File qnsort_bam = "~{outbam}"
    }

    String outbam = basename(bam, ".bam") + "picard-queryname-sorted.bam"

    String disk_type = if defined(num_ssds) then " LOCAL" else " SSD"

    Float N = size(bam, "GiB")
    Int scaleup_factor = if (N > 100) then 6 else 4
    Int persistend_disk_size = 20 + ceil(scaleup_factor * N)

    Int disk_size = if defined(num_ssds) then 375*select_first([num_ssds]) else persistend_disk_size

    command <<<
        set -eux

        # higher memory, also lower # of reads in memory given ~100 longer reads (1.5E4 bp vs 1.5E2 bp)
        gatk SortSam \
            --java-options "-Xmx28G -Xms24G" \
            -use_jdk_deflater -use_jdk_inflater \
            --MAX_RECORDS_IN_RAM 5000 \
            -I ~{bam} \
            -O ~{outbam} \
            -SO queryname
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          6,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.4.0.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + disk_type
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

#################################################
# intensive transformations -- split
#################################################

task SplitNameSortedUbam {
    meta {
        desciption: "Split a read-name-sorted unaligned BAM into chunks."
    }
    parameter_meta {
        read_cnt: "number of reads in the uBAM; providing this will reduce run time."
        n_reads: "desired number of reads per split; mutually exclusive with n_files"
        n_files: "desired number of split files; mutually exclusive with n_reads"
        uBAM: { localization_optional: true }
    }
    input {
        File uBAM
        Int? read_cnt
        Int? n_reads
        Int? n_files

        RuntimeAttr? runtime_attr_override
    }
    output {
        Array[File] split = glob("split_outputs/*.bam")
    }

    Boolean fail = defined(n_reads) == defined(n_files)  # mutually exclusive

    Int X = select_first([n_reads, n_files])
    String split_arg = if defined(n_reads) then "--SPLIT_TO_N_READS ~{X}" else "--SPLIT_TO_N_FILES ~{X}"
    String helper_arg = if (defined(read_cnt)) then "--TOTAL_READS_IN_INPUT ~{read_cnt}" else " "

    String base = basename(uBAM, ".bam")
    String local_bam = "/mnt/disks/cromwell_root/~{base}.bam"

    command <<<
        set -eux

        if ~{fail}; then echo "one and only one of [n_reads, n_files] must be specified" && exit 1; fi

        # prep
        time gcloud storage cp ~{uBAM} ~{local_bam}
        mkdir -p split_outputs

        # higher memory, also lower # of reads in memory given ~100 longer reads (1.5E4 bp vs 1.5E2 bp)
        gatk SplitSamByNumberOfReads \
            --java-options "-Xmx28G -Xms24G" \
            -use_jdk_deflater -use_jdk_inflater \
            --MAX_RECORDS_IN_RAM 5000 \
            -I ~{local_bam} \
            -O split_outputs \
            ~{split_arg} \
            ~{helper_arg}
    >>>
    #########################
    Int disk_size = 20 + ceil(3 * size(uBAM, "GiB"))

    RuntimeAttr default_attr = object {
        cpu_cores:          6,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.4.0.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SplitByRG {
    meta {
        description: "Split a BAM file that was aggregated, for the same sample, into pieces by read group."
    }
    input {
        File bam

        String out_prefix

        Int? num_ssds

        Boolean retain_rgless_records = false
        Boolean sort_and_index = false
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: "BAM to be split"
        out_prefix: "prefix for output bam and bai file names"
        sort_and_index: "if the user wants to (pos-)sort and index the resulting BAMs; this indicates the input BAM is mapped"

        split_bam: "the resuling BAMs, each having reads only in a single read group"
        split_bai: "the accompanying BAIs, if possible and explicit requested"
    }

    Int disk_size = if defined(num_ssds) then 375*select_first([num_ssds]) else 1+3*ceil(size([bam], "GB"))

    Array[String] extra_args = if (retain_rgless_records) then ["-u", "~{out_prefix}_noRG.bam"] else [""]
    command <<<
        set -eux

        samtools view -H ~{bam} | grep "^@RG" > "read_groups_header.txt"
        cat "read_groups_header.txt" | tr '\t' '\n' | grep "^ID:"  | awk -F ':' '{print $2}' > "RG_ids.txt"

        samtools split -@3 \
            -f "~{out_prefix}_%#.bam" \
            ~{sep=" " extra_args} \
            ~{bam}
        if ~{sort_and_index} ;
        then
            # cleanup space for the sorting
            rm ~{bam}
            for split_bam in "~{out_prefix}_"*.bam;
            do
                mv "${split_bam}" temp.bam
                samtools sort \
                    --write-index \
                    -o "${split_bam}##idx##${split_bam}.bai" \
                    temp.bam
            done
        fi
    >>>

    output {
        File read_group_header = "read_groups_header.txt"
        Array[String] rg_ids   = read_lines("RG_ids.txt")
        Array[File]  split_bam = glob("*.bam")
        Array[File?] split_bai = glob("*.bai")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ShardAlignedBam {
    meta {
        desciption: "Split an WGS BAM based on a provided scatter scheme."
    }
    parameter_meta {
        aligned_bam: {
            localization_optional: true,
            description: "input BAM file (must be coordinate sorted)."
        }
        aligned_bai: "input BAM index file"

        scatter_scheme: "A txt file holding how to scatter the WGS bam. Example (this example size-balance among the shards): ...\nchr5,chr19\nchr6,chrY,chrM\n..."

        parallel_subset_jobs: "an optimization; increasing this will lead to renting more powerfull VMs from GCP, though with shorter wall-clock time."
    }
    input {
        File  aligned_bam
        File? aligned_bai
        File scatter_scheme

        Int parallel_subset_jobs = 7  # empirical

        RuntimeAttr? runtime_attr_override
    }
    output {
        File unmapped_reads     = "~{base}.unmapped-reads.bam"

        Array[File] split_bams  = glob("~{base}.shard-*.bam")
    }

    Int disk_size = 3 * ceil(size(aligned_bam, "GB"))

    String base = basename(aligned_bam, ".bam")

    String local_bam = "/mnt/disks/cromwell_root/~{base}.bam"
    String local_bai = "/mnt/disks/cromwell_root/~{base}.bam.bai"

    Int vm_cores = parallel_subset_jobs * 2 + 2
    Int vm_memory = vm_cores * 4

    command <<<
        set -eux

        # here we use an optimization, that is, in stead of relying on the slow Cromwell localization,
        # we explicity localize the bam in the with gcloud storage cp
        time gcloud storage cp ~{aligned_bam} ~{local_bam}

        echo "==========================================================="
        echo "verify input bam is sorted by coordinate"
        samtools view -H ~{local_bam} | grep "@HD" > hd.line
        if ! grep -qF "SO:coordinate" hd.line;
        then
            echo "BAM must be coordinate sorted!" && echo && cat hd.line && exit 1
        fi
        echo "==========================================================="
        echo "index if bai not provided"
        if ~{defined(aligned_bai)}; then
            mv ~{aligned_bai} ~{local_bai}
        else
            time samtools index -@3 "~{local_bam}"
        fi
        echo "==========================================================="
        echo "######################################"
        echo "handle unmapped reads, if any, here"
        samtools view -@3\
            -f4 \
            -o "~{base}.unmapped-reads.bam" \
            "~{local_bam}" &
        echo "######################################"
        echo "first pad the provided sharding scheme with the uncovered contigs in the bam header"
        samtools view -H ~{local_bam} | grep "^@SQ" | awk -F '\t' '{print $2}' | awk -F ':' '{print $2}' > contigs.in.header.txt
        comm -13 \
            <(tr ',' '\n' < ~{scatter_scheme} | sort) \
            <(sort contigs.in.header.txt) \
            | tr '\n' ',' \
        > uncovered.scatter_scheme.txt
        cat uncovered.scatter_scheme.txt
        cat uncovered.scatter_scheme.txt >> ~{scatter_scheme}
        cat ~{scatter_scheme}
        echo "######################################"
        echo "now split according to the sharding scheme provided"
        job_cnt=0 # assume few unmapped reads, so don't count that
        idx=1
        while IFS= read -r one_shard
        do
            XX=$(echo "${one_shard}" | tr ',' ' ')
            read -ra YY <<< "$XX"
            samtools view -@1 \
                -o "~{base}.shard-${idx}.bam" \
                "~{local_bam}" \
                "${YY[@]}" &
            idx=$(( idx + 1 ))
            job_cnt=$(( job_cnt + 1 ))
            # let's not shoot ourselves
            if [[ ${job_cnt} == ~{parallel_subset_jobs} ]]; then wait; job_cnt=0; fi
        done < ~{scatter_scheme}
        wait
        echo "==========================================================="
        echo "DONE!"
        ls
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          vm_cores,
        mem_gb:             vm_memory,
        disk_gb:            disk_size,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

#################################################
# intensive transformations -- merge
#################################################

task MergeBamsWithSamtools {
    meta {
        description : "Merge several input BAMs into a single BAM."
        warn: "assumes input BAMs are coordinate sorted"
    }

    parameter_meta {
        bams: {localization_optional: true}
        out_prefix: "result file will be named <out_prefix>.bam"
    }

    input {
        Array[File] bams
        String out_prefix = "out"

        String disk_type = "LOCAL"

        RuntimeAttr? runtime_attr_override
    }

    output {
        File merged_bam = "~{out_prefix}.bam"
        File merged_bai = "~{out_prefix}.bam.bai"
    }

    command <<<
        set -euxo pipefail

        mkdir -p bams_dir
        time \
        gcloud storage cp \
            ~{sep=' ' bams} \
            /mnt/disks/cromwell_root/bams_dir/
        ls bams_dir

        cd bams_dir && ls ./*.bam > bams.list
        time \
        samtools merge \
            -p -c --no-PG \
            -@ 3 \
            --write-index \
            -o "~{out_prefix}.bam##idx##~{out_prefix}.bam.bai" \
            -b bams.list
        mv ~{out_prefix}.bam \
           ~{out_prefix}.bam.bai \
        /mnt/disks/cromwell_root/
    >>>
    #########################
    Int local_ssd_sz = if size(bams, "GiB") > 150 then 750 else 375
    Int pd_sz = 10 + 3*ceil(size(bams, "GiB"))
    Int disk_size = if "LOCAL" == disk_type then local_ssd_sz else pd_sz

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        preemptible_tries:  if "LOCAL" == disk_type then 1 else 0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MergeBamsQuerynameSortedWithPicard {
    meta {
        desciption: "Merge list of bams that were queryname sorted with Picard"
    }
    parameter_meta {
        qns_bams: {
            desciption: "queryname-sorted, preferrably by Picard, bams to be merged",
            localization_optional: true
        }
        base_names: "basenames of all files, INCLUDING the '.bam' extention."
        out_prefix: "result file will be named <out_prefix>.bam"
        num_ssds: "if provided, will use LOCAL SSDs for faster speed at higher cost"
    }
    input {
        Array[File] qns_bams
        Array[String] base_names
        String out_prefix

        Int? num_ssds
        RuntimeAttr? runtime_attr_override
    }
    output {
        File res = "~{out_prefix}.bam"
    }

    Float N = ceil(size(qns_bams, "GB"))
    Int scaleup_factor = if (N > 100) then 6 else 4
    Int persistend_disk_size = 20 + ceil(scaleup_factor * N)

    Int disk_size = if defined(num_ssds) then 375*select_first([num_ssds]) else persistend_disk_size
    String disk_type = if defined(num_ssds) then " LOCAL" else " SSD"

    command <<<
        set -eux

        mkdir -p bams_dir
        gcloud storage cp \
            ~{sep=' ' qns_bams} \
            /mnt/disks/cromwell_root/bams_dir/
        ls bams_dir

        # higher memory, also lower # of reads in memory given ~100 longer reads (1.5E4 bp vs 1.5E2 bp)
        cd bams_dir
        gatk MergeSamFiles \
            --java-options "-Xmx28G -Xms24G" \
            --USE_THREADING \
            -use_jdk_deflater -use_jdk_inflater \
            --SORT_ORDER queryname \
            -I ~{sep=" -I " base_names} \
            -O "/mnt/disks/cromwell_root/~{out_prefix}.bam"
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          6,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-custom-gatk:4.4.0.0-samtools1.18"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + disk_type
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
