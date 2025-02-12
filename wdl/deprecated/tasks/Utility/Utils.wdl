version 1.0

import "../../../structs/Structs.wdl"

task SplitDelimitedString {

    meta {
        description: "Split a delimited string into an array"
    }

    parameter_meta {
        s: "The string to split"
        separate: "The delimiter to split on"
    }

    input {
        String s
        String separate
    }

    command <<<
        set -eux

        echo ~{s} | tr ~{separate} '\n' > result.txt
    >>>

    output {
        Array[String] arr = read_lines("result.txt")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
