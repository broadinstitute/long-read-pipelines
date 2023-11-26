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

        tmp_tsv=~{write_map(input_map)}
        awk -F '\t' '{print $1}' ${tmp_tsv} > keys.txt
        awk -F '\t' '{print $2}' ${tmp_tsv} > values.txt
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
