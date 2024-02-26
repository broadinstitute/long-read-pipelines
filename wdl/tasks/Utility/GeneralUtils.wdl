version 1.0

# todo: move all existing, simiar util tasks here
# Hosting utils that just needs basic Linux shell programs

task TarGZFiles {
    meta {
        description:
        "Zip up a list of files to a tar.gz file."
    }

    parameter_meta {
        files: "List of files to zip up."
        name: "Name of the tar.gz file."
    }

    input {
        Array[File] files
        String name
    }

    command <<<
        set -eux
        mkdir -p save/
        for ff in ~{sep=' ' files}; do cp "${ff}" save/ ; done
        tar -cvzf ~{name}.tar.gz -C save/ .
    >>>

    output {
        File you_got_it = "~{name}.tar.gz"
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task GetTodayDate {
    meta {
        desciption: "Generates a YYYY-MM-DD date of today (when this task is called). UTC."
        volatile: true
    }
    command {
        date '+%Y-%m-%d'
    }

    output {
        String yyyy_mm_dd = read_string(stdout())
    }
    runtime {
        disks: "local-disk 10 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task GetBetterRange {
    meta {
        desciption:
        "Given the start, step-size, and end value, return an array of numbers"
        note:
        "Basically I'm a linux `seq` command; I exist since WDL stdlib range function has limits."
    }

    input {
        Int start
        Int end
        Int step
    }

    output {
        Array[Int] your_range = read_lines("result.txt")
    }
    #############################################################
    command <<<
        set -euxo pipefail
        seq ~{start} ~{step} ~{end} | tee > "result.txt"
    >>>

    #############################################################
    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
        preemptible: 1
        maxRetries: 1
    }
}

task CollapseArrayOfStrings {
    meta {
        desciption: "For collapsing an array of strings using space."
        note: "When the next version (> 1.0) of WDL is supported on Terra, use the official solution."
    }
    input {
        Array[String] input_array
        String joiner
    }
    output {
        String collapsed = read_string("result.txt")
    }

    command <<<
        set -euxo pipefail

        n=$(echo ~{joiner} | wc -c | awk '{print $1}')
        if [[ ! "${n}" -eq 1 ]]; then echo "cannot collapse with multi-char joiner" && exit 1; fi

        tr '\n' "~{joiner}" < ~{write_lines(input_array)} \
        > result.txt
    >>>
    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task CoerceMapToArrayOfPairs {
    meta {
        desciption: "To coerce a WDL Map into Array[Pair], since Cromwell doesn't support it. Mostly used for iterating a map."
    }

    input {
        Map[String, String] input_map
    }

    command <<<
        set -eux

        two_col_tsv=~{write_map(input_map)}
        cat "${two_col_tsv}"
        wc -l "${two_col_tsv}"
        # because some Cromwell versions' stdlib function write_map() doesn't have new line at end of file, so we add it explicitly
        if [[ $(tail -c1 "${two_col_tsv}" | wc -l) -eq 0 ]]; then
            sed -i -e '$a\' "${two_col_tsv}"
        fi
        # '
        wc -l "${two_col_tsv}"
        awk -F '\t' '{print $1}' ${two_col_tsv} > keys.txt
        awk -F '\t' '{print $2}' ${two_col_tsv} > values.txt
    >>>

    output {
        Array[Pair[String, String]] output_pairs = zip(read_lines("keys.txt"), read_lines("values.txt"))
    }

    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task CoerceArrayOfPairsToMap {
    meta {
        description:
        "To coerce an Array of Pair's to a Map; use only when you're sure the 'key' array is unique."
    }
    input {
        Array[String] keys
        Array[String] values
    }

    command <<<
        set -eux

        paste ~{write_lines(keys)} ~{write_lines(values)} > "result.tsv"
    >>>

    output {
        Map[String, String] output_map = read_map("result.tsv")
    }
    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task MergeMaps {
    meta {
        desciption:
        "For merging two maps into one."
        note:
        "User is responsible for ensuring uniqueness of keys."
    }
    input {
        Map[String, String] one
        Map[String, String] two
    }
    output {
        Map[String, String] merged = read_map("merged.tsv")
    }

    command <<<
        set -euxo pipefail

        ####### one
        two_col_tsv=~{write_map(one)}
        cat "${two_col_tsv}"
        wc -l "${two_col_tsv}"
        # because some Cromwell versions' stdlib function write_map() doesn't have new line at end of file, so we add it explicitly
        if [[ $(tail -c1 "${two_col_tsv}" | wc -l) -eq 0 ]]; then
            sed -i -e '$a\' "${two_col_tsv}"
        fi
        # '
        wc -l "${two_col_tsv}"
        mv "${two_col_tsv}" "one.tsv"
        ####### two
        two_col_tsv=~{write_map(two)}
        cat "${two_col_tsv}"
        wc -l "${two_col_tsv}"
        # because some Cromwell versions' stdlib function write_map() doesn't have new line at end of file, so we add it explicitly
        if [[ $(tail -c1 "${two_col_tsv}" | wc -l) -eq 0 ]]; then
            sed -i -e '$a\' "${two_col_tsv}"
        fi
        # '
        wc -l "${two_col_tsv}"
        mv "${two_col_tsv}" "two.tsv"
        ####### merge
        cat "one.tsv" "two.tsv" > "merged.tsv"
        cat "merged.tsv"
    >>>

    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task MapToTsv {
    meta {
        desciption:
        "For fixing an issue of some Cromwell servers where write_map misses the last line"
    }
    input {
        Map[String, String] m
    }
    output {
        File tsv = "two_col.tsv"
    }
    command <<<
        two_col_tsv=~{write_map(m)}
        cat "${two_col_tsv}"
        wc -l "${two_col_tsv}"
        # because some Cromwell versions' stdlib function write_map() doesn't have new line at end of file, so we add it explicitly
        if [[ $(tail -c1 "${two_col_tsv}" | wc -l) -eq 0 ]]; then
            sed -i -e '$a\' "${two_col_tsv}"
        fi
        # '
        wc -l "${two_col_tsv}"
        mv "${two_col_tsv}" "two_col.tsv"
    >>>
    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task ConcatenateFiles {
    meta {
        desciption:
        "For concatinating files"
    }
    input {
        Array[File]+ af
        String out_name
    }
    output {
        File merged = "~{out_name}"
    }
    command <<<
        set -euxo pipefail

        cat ~{sep=' ' af} > "~{out_name}"
    >>>
    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

struct HelperStructForUnzip {
    Array[Pair[String, String]] contents
}
task Unzip {
    meta {
        desciption:
        "unzip an array of pairs to a pair of two arrays of string"
        note:
        "this task also serves as how to do jq operations for in-task unzip"
    }
    input {
        Array[Pair[String, String]] apss
    }
    output {
        Pair[Array[String], Array[String]] res = (read_lines("left.txt"), read_lines("right.txt"))
    }

    HelperStructForUnzip x = object {contents: apss}
    command <<<
    set -euxo pipefail

        mv ~{write_json(x)} tmp.json

        jq --raw-output '.contents[] | .left' tmp.json  > left.txt
        jq --raw-output '.contents[] | .right' tmp.json > right.txt

        cat left.txt
        cat right.txt
    >>>

    runtime {
        disks: "local-disk 10 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/jq:1.7.1"
    }
}

task SendEmailNotification {
    meta {
        desciption:
        "Send email, e.g. when certain events happen."
        warn:
        "This may exhaust your SendGrid API quota, use caution or pay up."
    }
    parameter_meta {
        sendgrid_api_key_file: "A JSON file holding the account key (guard it carefully)"
        sender_name: "Name of the sender of the email"
        sender_email: "Email address registered on SendGrid used to send email with."

        receiver_names_and_addresses: "intended receivers (don't spam)"

        email_subject: "The subject/title/topic of the email"
        email_body: "The plain-text contents of the email"

        html: "a single HTML document to include in the email"

        txt_attachment_names_and_files: "Names and Files of .txt files to attach to the email."
        tsv_attachment_names_and_files: "Names and Files of .tsv files to attach to the email."
        pdf_attachment_names_and_files: "Names and Files of .pdf files to attach to the email."
    }
    input {
        File sendgrid_api_key_file
        String sender_name
        String sender_email

        Array[Pair[String, String]] receiver_names_and_addresses

        String email_subject
        String email_body

        File? html

        Array[Pair[String, File]]? txt_attachment_names_and_files
        Array[Pair[String, File]]? tsv_attachment_names_and_files
        Array[Pair[String, File]]? pdf_attachment_names_and_files
    }

    HelperStructForUnzip receivers = object {contents: receiver_names_and_addresses}

    Boolean has_txt_attach = defined(txt_attachment_names_and_files)
    HelperStructForUnzip txt_attach_obj = object {contents: select_first([txt_attachment_names_and_files, [('null', 'null')]])}

    Boolean has_tsv_attach = defined(tsv_attachment_names_and_files)
    HelperStructForUnzip tsv_attach_obj = object {contents: select_first([tsv_attachment_names_and_files, [('null', 'null')]])}

    Boolean has_pdf_attach = defined(pdf_attachment_names_and_files)
    HelperStructForUnzip pdf_attach_obj = object {contents: select_first([pdf_attachment_names_and_files, [('null', 'null')]])}

    command <<<
    set -euxo pipefail

        ##########################################################
        # some boiler-plate stuff to just do arg massaging
        ##########################################################
        ## receivers
        mv ~{write_json(receivers)} tmp.json

        jq --raw-output '.contents[] | .left' tmp.json  > receiver_names.txt
        jq --raw-output '.contents[] | .right' tmp.json > receiver_emails.txt
        rm tmp.json
        ## txt
        if ~{has_txt_attach}; then
            mv ~{write_json(txt_attach_obj)} tmp.json
            bash /opt/localize_files.sh \
                tmp.json \
                $(pwd) \
                txt.attach.tsv
            rm tmp.json
        fi
        ## tsv
        if ~{has_tsv_attach}; then
            mv ~{write_json(tsv_attach_obj)} tmp.json
            bash /opt/localize_files.sh \
                tmp.json \
                $(pwd) \
                tsv.attach.tsv
            rm tmp.json
        fi
        ## pdf
        if ~{has_pdf_attach}; then
            mv ~{write_json(pdf_attach_obj)} tmp.json
            bash /opt/localize_files.sh \
                tmp.json \
                $(pwd) \
                pdf.attach.tsv
            rm tmp.json
        fi

        ##########################################################
        # kick off
        ##########################################################
        python3 /opt/send_email.py \
            --sendgrid_api_key ~{sendgrid_api_key_file} \
            --sender_name ~{sender_name} \
            --sender_email ~{sender_email} \
            --notification_receiver_names receiver_names.txt \
            --notification_receiver_emails receiver_emails.txt \
            --email_subject "~{email_subject}" \
            --email_body "~{email_body}" \
            ~{true='--txt_names_and_files' false=' ' has_txt_attach} \
            ~{true='txt.attach.tsv'        false=' ' has_txt_attach} \
            ~{true='--tsv_names_and_files' false=' ' has_tsv_attach} \
            ~{true='tsv.attach.tsv'        false=' ' has_tsv_attach} \
            ~{true='--pdf_names_and_files' false=' ' has_pdf_attach} \
            ~{true='pdf.attach.tsv'        false=' ' has_pdf_attach}

    >>>
    runtime {
        disks: "local-disk 10 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-wdl-email:0.0.1"
    }
}
