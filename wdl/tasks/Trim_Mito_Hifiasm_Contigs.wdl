version 1.0
import "Structs.wdl"

workflow Split_fasta {
    input {
        File fasta
    }
    call Self_Align {
        input:
        fasta = fasta
    }
    output{
        Array[split_fasta] output_fasta = Self_Align.split_fasta
    }
}


task Self_Align {
    input {
        File fasta
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        fasta: "Hifiasm Assembly Fasta File"

    }

    Int disk_size = 50 * ceil(size(fasta, "GB"))

    command <<<
        set -euxo pipefail

        while read line ; do
            if [ ${line:0:1} == ">" ]; then
                filename=$(echo "$line" | cut -d ":" -f1 | tr -d ">")
                touch ./"$filename".fasta
                echo "$line" >> ./"${filename}".fasta
            else
                echo "$line" >> ./"${filename}".fasta
            fi
        done < ${fasta}
    >>>

    output {
        Array[File] split_fasta = "split_fasta.fa"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-canu:0.2.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}