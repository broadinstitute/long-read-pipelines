version 1.0

import "../../../structs/Structs.wdl"

workflow PfalciparumHrp2Hrp3Status {

    meta {
        description: "Determine if HRP2 and HRP3 are deleted in a sample"
    }
    parameter_meta {
        bam: "BAM file"
        bai: "BAM index file"
    }

    input {
        File bam
        File bai
    }

    call IsLocusDeleted as HRP2Status {
        input:
            bam = bam,
            bai = bai,
            chr = "Pf3D7_08_v3",
            start = 1373212,
            stop = 1376988
    }

    call IsLocusDeleted as HRP3Status {
        input:
            bam = bam,
            bai = bai,
            chr = "Pf3D7_13_v3",
            start = 2840236,
            stop = 2842840
    }

    # Finalize data

    output {
        String hrp2 = HRP2Status.status
        String hrp3 = HRP3Status.status
    }
}

task IsLocusDeleted {
    input {
        File bam
        File bai

        String chr
        Int start
        Int stop

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(bam, "GB"))  + 2*ceil(size(bai, "GB"))

    command <<<
        set -euxo pipefail

        samtools view -bhX ~{bam} ~{bai} ~{chr} > chr.bam
        samtools index chr.bam

        mosdepth -t 4 -b <(echo -e "~{chr}\t~{start}\t~{stop}") -x -Q 1 out chr.bam

        grep total out.mosdepth.summary.txt | \
            paste - - | \
            awk '{ if ($10 > 0.05*$4) print "+"; else print "-" }' \
            > status.txt
    >>>

    output {
        String status = read_string("status.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-mosdepth:0.3.1"
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
