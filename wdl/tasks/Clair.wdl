version 1.0

#################################################
# This pipeline calls small variants using Clair.
#################################################

import "Structs.wdl"

task Clair {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fai

        String model_class
        String chr

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: "input BAM from which to call variants"
        bai: "index accompanying the BAM"

        ref_fasta: "reference to which the BAM was aligned to"
        ref_fai:   "index accompanying the reference"

        model_class: "class of model to be applied"
        chr: "chromsome on which to call variants"
    }

    Int num_cpus = 4
    Int disk_size = ceil(size(bam, "GB")) + 50
    String prefix = basename(bam, ".bam")

    command <<<
        set -euxo pipefail

        SAMPLE_NAME=$(samtools view -H ~{bam} | grep -m1 '@RG' | sed 's/\t/\n/g' | grep 'SM' | awk -F":" '{ print $2 }')
        MODEL=$(echo ~{model_class} | tr 'A-Z' 'a-z')

        python3 $CLAIR callVarBam \
            --chkpnt_fn /opt/clair/models/$MODEL/model \
            --ref_fn ~{ref_fasta} \
            --bam_fn ~{bam} \
            --ctgName ~{chr} \
            --threads ~{num_cpus} \
            --sampleName $SAMPLE_NAME \
            --call_fn ~{prefix}.clair.~{chr}.vcf
    >>>

    output {
        File vcf = "~{prefix}.clair.~{chr}.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-clair:2.1.1"
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
