version 1.0

import "../../structs/Structs.wdl"

task FastP {
    input {
        File illumina_fq1
        File illumina_fq2

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        2,
        docker:             "quay.io/biocontainers/fastp:0.23.2--h5f740d0_3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int num_cpu = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    command <<<
        set -euxo pipefail

        mkdir output
        fastp --in1 ~{illumina_fq1} --in2 ~{illumina_fq2} \
            --out1 "output/~{basename(illumina_fq1)}" --out2 "output/~{basename(illumina_fq2)}" \
            --unpaired1 "output/unpaired.fq.gz" --unpaired2 "output/unpaired.fq.gz" \
            --html output/fastp.html --json output/fastp.json \
            --threads ~{num_cpu}
    >>>

    output {
        File processed_fq1 = "output/~{basename(illumina_fq1)}"
        File processed_fq2 = "output/~{basename(illumina_fq2)}"
        File unpaired_fq = "output/unpaired.fq.gz"
        File fastp_report = "output/fastp.html"
        File fastp_json = "output/fastp.json"
    }
}