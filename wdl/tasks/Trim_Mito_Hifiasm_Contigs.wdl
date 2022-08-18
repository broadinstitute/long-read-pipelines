version 1.0


workflow Trim_Contigs {
    input {
        File assembly_fasta
    }

    call Filter_Contigs {
        input:
        assembly_fasta = assembly_fasta
    }

#    call Self_Align {
#        input:
#        assembly_fasta = assembly_fasta
#    }

    output{
        #Array[File] output_fasta = Self_Align.split_fasta
        File selected_contigs = Filter_Contigs.filtered_contigs
    }
}


task Filter_Contigs {
    input {
        File assembly_fasta
    }

    parameter_meta {
        assembly_fasta: "Hifiasm Assembly Fasta File"
    }

    command <<<
        set -euxo pipefail
        awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >16000 && length($2) {print ">"$0}' ~{assembly_fasta} > filtered_contigs.fasta

    >>>

    output {
        File filtered_contigs = "filtered_contigs.fasta"
    }

    #########################
    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}


task Self_Align {
    input {
        File assembly_fasta
    }

    parameter_meta {
        assembly_fasta: "Hifiasm Assembly Fasta File"

    }



    command <<<
        set -euxo pipefail

        while read line ; do
            if [ ${line:0:1} == ">" ]; then
                filename=$(echo "$line" | cut -d ":" -f1 | tr -d ">")
                touch "$filename".fasta
                echo "$line" >> "split.${filename}".fasta
            else
                echo "$line" >> "split.${filename}".fasta
            fi
        done < ~{assembly_fasta}
    >>>

    output {
        Array[File] split_fasta = glob("split.*.fasta")
    }


    #########################
    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}