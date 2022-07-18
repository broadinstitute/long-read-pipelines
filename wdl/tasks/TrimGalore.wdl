version 1.0

import "Structs.wdl"

task TrimGalore {
    input {
        File reads_fq1
        File? reads_fq2

        Int? min_length
        RuntimeAttr? runtime_attr_override
    }

    Int num_cores = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    Array[File] reads = select_all([reads_fq1, reads_fq2])
    Boolean paired = length(reads) > 1
    Int length = select_first([min_length, 31])

    command <<<
        mkdir output
        trim_galore -t ~{num_cores} --length ~{length} ~{true='--paired' false='' paired} ~{sep=' ' reads} -o output
    >>>

    output {
        File trimmed_fq1 = "output/~{basename(reads_fq1)}"
        File trimming_report1 = "output/~{basename(reads_fq1)}_trimming_report.txt"
        File? trimmed_fq2 = "output/~{basename(reads_fq2)}"
        File? trimming_report2 = "output/~{basename(reads_fq2)}_trimming_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            ceil((size(reads_fq1, "G") + size(reads_fq2, "G")) * 2 + 10.0),
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0"
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
