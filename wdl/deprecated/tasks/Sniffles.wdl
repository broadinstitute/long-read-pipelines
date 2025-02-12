version 1.0

import "../../structs/Structs.wdl"

# Given BAM, call SVs using Sniffles
task Sniffles {

    meta {
        description: "Call SVs using Sniffles-1"
    }

    parameter_meta {
        bam:              "input BAM from which to call SVs"
        bai:              "index accompanying the BAM"
        min_read_support: "[default-valued] minimum reads required to make a call"
        min_read_length:  "[default-valued] filter out reads below minimum read length"
        min_mq:           "[default-valued] minimum mapping quality to accept"
        chr:              "chr on which to call variants"
        prefix:           "prefix for output file"
        runtime_attr_override: "override default runtime attributes"
    }

    input {
        File bam
        File bai
        Int min_read_support = 2
        Int min_read_length = 1000
        Int min_mq = 20
        String? chr
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    Int disk_size = 2*ceil(size([bam, bai], "GB"))
    String fileoutput = if defined(chr) then "~{prefix}.~{chr}.sniffles.vcf" else "~{prefix}.sniffles.vcf"

    command <<<
        set -x

        sniffles -t ~{cpus} \
                 -m ~{bam} \
                 -v ~{fileoutput}\
                 -s ~{min_read_support} \
                 -r ~{min_read_length} \
                 -q ~{min_mq} \
                 --num_reads_report -1 \
                 --genotype

        touch ~{prefix}.~{fileoutput}

        cat ~{prefix}.~{fileoutput}| \
            grep -v -e '##fileDate' | \
            awk '{ if ($1 ~ "^#" || $7 == "PASS") print $0 }' \
            > ~{prefix}.~{fileoutput}
    >>>

    output {
        File vcf = "~{fileoutput}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             46,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sv:0.1.8"
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

