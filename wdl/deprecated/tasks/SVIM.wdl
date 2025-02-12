version 1.0

##########################################################################################
# This pipeline calls SVs on an input LR BAM using various known SV algorithms
# that are specifically designed to work with long read data.
# Each individual task/algo. is directly callable, if so desired.
##########################################################################################

import "../../structs/Structs.wdl"

# Given BAM, call SVs using SVIM
task SVIM {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: "input BAM from which to call SVs"
        bai: "index accompanying the BAM"

        ref_fasta: "reference to which the BAM was aligned to"
        ref_fasta_fai: "index accompanying the reference"

        prefix: "prefix for output"
    }

    Int disk_size = 2 * ceil(size([bam, bai, ref_fasta, ref_fasta_fai], "GB"))

    command <<<
        set -euxo pipefail

        SM=$(samtools view -H ~{bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g')

        svim alignment \
            --sample ${SM} \
            --insertion_sequences \
            --read_names ~{prefix}_svim_files \
            ~{bam} \
            ~{ref_fasta}

        grep -v -e '##fileDate' ~{prefix}_svim_files/variants.vcf > ~{prefix}.svim.vcf
    >>>

    output {
        File vcf = "~{prefix}.svim.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
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
