version 1.0

import "../../structs/Structs.wdl"

task CollapseIsoforms {
    input {
        File bam
        File ref_fasta

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(2*size(ref_fasta, "GB") + 10*size(bam, "GB"))

    command <<<
        set -euxo pipefail

        source activate lr-tama

        samtools view -h -F 4 ~{bam} > tmp.sam
        python /tama/tama_collapse.py -s tmp.sam -f ~{ref_fasta} -p ~{prefix} -x capped
    >>>

    output {
        File collapsed_read = "~{prefix}_read.txt"
        File collapsed_variants = "~{prefix}_variants.txt"
        File collapsed_polya = "~{prefix}_polya.txt"
        File collapsed_strand_check = "~{prefix}_strand_check.txt"
        File collapsed_collapsed_varcov = "~{prefix}_varcov.txt"
        File collapsed_trans_report = "~{prefix}_trans_report.txt"
        File collapsed_collapsed = "~{prefix}.bed"
        File collapsed_local_density_error = "~{prefix}_local_density_error.txt"
        File collapsed_trans_read = "~{prefix}_trans_read.bed"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-tama:0.1.1"
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
