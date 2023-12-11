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
