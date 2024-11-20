version 1.0

import "../../structs/Structs.wdl"

workflow PhaseReads {
    input {
        File reads
        String PL
        Boolean? reads_are_corrected
    }

    call ConvertToFastq {
        input:
            reads = reads
    }

    call Minimap2PAF {
        input:
            reads_fq = ConvertToFastq.reads_fq,
            PL = PL
    }
}

task ConvertToFastq {
    input {
        File reads

        RuntimeAttr? runtime_attr_override
    }

    String fqname = basename(reads, ".bam") + ".fastq.gz"
    Int disk_size = 4*ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        samtools fastq ~{reads} | gzip -1 > ~{fqname}
    >>>

    output {
        File reads_fq = fqname
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
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

task Minimap2PAF {
    input {
        File reads_fq
        String PL
        Boolean? reads_are_corrected

        RuntimeAttr? runtime_attr_override
    }

    Boolean correct = select_first([reads_are_corrected, false])
    String map_preset = if (PL == "ONT") then "ava-ont" else "ava-pb"

    Int cpus = 4
    Int disk_size = 4*ceil(size(reads_fq, "GB"))
    String paf_name = sub(reads_fq, ".fastq.gz", ".paf.gz")

    command <<<
        set -euxo pipefail

        minimap2 -cx ~{map_preset} -t ~{cpus} ~{reads_fq} ~{reads_fq} | gzip -1 > ~{paf_name}
    >>>

    output {
        File paf = paf_name
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
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
