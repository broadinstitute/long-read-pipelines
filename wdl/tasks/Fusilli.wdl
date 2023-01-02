version 1.0

import "Structs.wdl"

task FusilliAssemble {
    input {
        File illumina_fq1
        File illumina_fq2

        Array[File] references

        File? previous_run

        Int k = 47

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + (3 * ceil(size(illumina_fq1) + size(illumina_fq2))) + ceil((length(references) * 25) / 1024)

    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us-east1-docker.pkg.dev/broad-dsp-lrma/fusilli/fusilli:devel"
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

    command <<<
        # Activate fusilli conda env
        source /usr/local/bin/_activate_current_env.sh
        set -euxo pipefail

        prev_run="~{default="" previous_run}"
        if [[ "${prev_run}" != "" ]]; then
            mkdir -p output
            tar -xzf "${prev_run}" -C output
        fi

        fusilli assemble -t ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            -k ~{k} \
            -R ~{write_lines(references)} \
            -c -o output/ ~{illumina_fq1} ~{illumina_fq2}

        tar -czf -C output/ fusilli_run.tar.gz .
    >>>

    output {
        File contigs = "output/contigs.fasta"
        File fusilli_output_tar = "fusilli_run.tar.gz"
    }
}


task FinalizeFusilliRun {
    input {
        File fusilli_output_tar
        String gcs_output_file

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        gsutil cp ~{fusilli_output_tar} ~{gcs_output_file}
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            20,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
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
