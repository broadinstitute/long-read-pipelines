version 1.0


import "tasks/Structs.wdl"

workflow Sniffles2 {
    input {
        Array[File] sampleBAMs
        Array[File] bai
        Int min_read_support
        Int min_alignment_length
        Int min_mq
        String prefix}

    scatter (bam in sampleBAMs) {
        call sample_sv {
                input:
                    bam = bam,
                    bai = bai,
                    min_read_support = min_read_support,
                    min_alignment_length = min_alignment_length,
                    min_mq = min_mq,
                    prefix = prefix
                }
         }

    call merge_call { input: snfs = sample_sv.snf}


    output {
         File multisample_vcf = merge_call.vcf
        }
}



# Given BAM, call SVs using Sniffles
task sample_sv {

    input {
        File bam
        Array[File] bai
        Int min_read_support #= 2
        Int min_alignment_length #= 1000
        Int min_mq #= 20
        String? chr
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:              "input BAM from which to call SVs"
        bai:              "index accompanying the BAM"
        min_read_support: "[default-valued] minimum reads required to make a call"
        min_alignment_length:  "[default-valued] reads shorter than the minimum alignment length will be ignored"
        min_mq:           "[default-valued] minimum mapping quality to accept"
        chr:              "chr on which to call variants"
        prefix:           "prefix for output"
    }

    Int cpus = 8
    Int disk_size = 2*ceil(size([bam, bai], "GB"))
    String fileoutput = if defined(chr) then "~{prefix}.~{chr}.sniffles.snf" else "~{prefix}.sniffles.snf"

    command <<<
        set -x

        sniffles -t ~{cpus} \
                 -i ~{bam} \
                 --minsupport ~{min_read_support} \
                 --min-alignment-length ~{min_alignment_length} \
                 --mapq ~{min_mq} \
                 --snf ~{fileoutput}
        tree
        touch ~{prefix}.~{fileoutput}

    >>>

    output {
        File snf = "~{fileoutput}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             46,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sniffles2:2.0.6"
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



task merge_call {
    input {
        Array[File] snfs
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
    snfs: "snf files"
}

    command <<<
    set -x
    sniffles --input ~{sep=" " snfs} \
    --vcf multisample.vcf
    tree

    >>>

    output {
        File vcf = "multisample.vcf" }

    Int cpus = 8
    Int disk_size = 3*ceil(size(snfs, "GB"))
#    Int disk_size = 60
                                                                                                                                                                                                                                                                                                                                                                                                                                   #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             46,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sniffles2:2.0.6"
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