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
