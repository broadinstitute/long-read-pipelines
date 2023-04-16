version 1.0

import "../../structs/Structs.wdl"

task Unicycler {
    input {
        String sample_name

        File illumina_fq1
        File illumina_fq2
        File? illumina_unpaired
        File? long_reads

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        2,
        docker:             "quay.io/biocontainers/unicycler:0.5.0--py310hc8f18ef_2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int num_cpu = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    command <<<
        set -euxo pipefail

        unicycler --threads ~{num_cpu} \
            -1 ~{illumina_fq1} -2 ~{illumina_fq2} \
            ~{'-s ' + illumina_unpaired} \
            ~{'-l ' + long_reads} \
            -o output

        cp "output/assembly.fasta" "~{sample_name}.fasta"
        cp "output/assembly.gfa" "~{sample_name}.gfa"
    >>>

    output {
        File assembly_fasta = "~{sample_name}.fasta"
        File assembly_gfa = "~{sample_name}.gfa"
        File unicycler_log = "output/unicycler.log"
    }

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