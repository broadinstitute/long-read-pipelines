version 1.0

import "Structs.wdl"

task BWAMem2Align {
    input {
        File ref
        Array[File]+ reads
        String output_prefix

        File? ref_0123 = ref + ".0123"
        File? ref_amb = ref + ".amb"
        File? ref_ann = ref + ".ann"
        File? ref_bwt = ref + ".bwt.2bit.64"
        File? ref_pac = ref + ".pac"

        Boolean? keep_unaligned = false

        RuntimeAttr? runtime_attr_override
    }

    output {
        File aligned_bam = "~{output_prefix}.bam"
        File aligned_bai = "~{output_prefix}.bam.bai"
    }

    command <<<
        set -euxo pipefail

        if [[ ! -f "~{ref_0123}" || ! -f "~{ref_amb}" || -f "~{ref_ann}" || ! -f "~{ref_bwt}" || ! -f "~{ref_pac}" ]]; then
            bwa-mem2 index "~{ref}"
        fi

        if [[ "~{true="1" false="0" keep_unaligned}" == "0" ]]; then
            bwa-mem2 mem -t ~{runtime_attr.cpu_cores} "~{ref}" ~{sep=' ' reads} \
                | samtools view -h -F 4 \
                | samtools sort -@ ~{runtime_attr.cpu_cores} -O bam -o "~{output_prefix}.bam"
        else
            bwa-mem2 mem -t ~{runtime_attr.cpu_cores} "~{ref}" ~{sep=' ' reads} \
                | samtools sort -@ ~{runtime_attr.cpu_cores} -O bam -o "~{output_prefix}.bam"
        fi

        samtools index "~{output_prefix}.bam"
    >>>


    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            20,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-bwa-mem2:2.2.1"
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

