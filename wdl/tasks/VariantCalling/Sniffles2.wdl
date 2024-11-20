version 1.0


import "../../structs/Structs.wdl"


workflow Sniffles2 {

    meta {
        description: "This workflow calls SV candidates using Sniffles2 population mode."
    }

    parameter_meta {
        # input
        sampleBAMs:      "GCS paths to aligned BAM files from multiple samples"
        sampleBAIs:       "GCS paths to aligned BAM files indices from multiple samples"
        sampleIDs:       "matching sample IDs of the BAMs"
        minsvlen:        "Minimum SV length in bp"
        prefix:          "prefix for output files"
        # output
        single_snf:      "[OUTPUT] .snf output containing SV candidates from a single sample"
        multisample_vcf: "[OUTPUT] Multi-sample vcf output"
    }

    input {
        Array[File] sampleBAMs
        Array[File] sampleBAIs
        Array[String] sampleIDs
        String prefix
        Int minsvlen = 50
    }

    scatter (i in range(length(sampleBAMs))) {
        call SampleSV {
            input:
                bam = sampleBAMs[i],
                bai = sampleBAIs[i],
                minsvlen = minsvlen,
                prefix = prefix,
                sample_id = sampleIDs[i]
        }
    }

    call MergeCall {
        input:
            snfs = SampleSV.snf,
            prefix = prefix
    }


    output {
         Array[File] single_snf = SampleSV.snf
         File multisample_vcf = MergeCall.vcf
    }
}



task SampleSV {

    meta {
        description: "This task calls SV candidates from a single sample."
    }

    input {
        File bam
        File bai
        Int minsvlen = 50
        String sample_id
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:              "input BAM from which to call SVs"
        bai:              "index accompanying the BAM"
        minsvlen:         "minimum SV length in bp. Default 50"
        sample_id:        "Sample ID"
        prefix:           "prefix for output"
    }

    Int cpus = 8
    Int disk_size = 2*ceil(size([bam, bai], "GB"))
    String snf_output = "~{prefix}.sniffles.snf"
    String vcf_output = "~{prefix}.sniffles.vcf"

    command <<<
        set -eux

        sniffles -t ~{cpus} \
                 -i ~{bam} \
                 --minsvlen ~{minsvlen} \
                 --sample-id ~{sample_id} \
                 --vcf ~{vcf_output} \
                 --snf ~{snf_output}
        tree
    >>>

    output {
        File snf = "~{snf_output}"
        File vcf = "~{vcf_output}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             46,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
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


task MergeCall {

    meta {
        description: "This tasks performs joined-calling from multiple .snf files and produces a single .vcf"
    }

    input {
        Array[File] snfs
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        snfs: ".snf files"
    }

    command <<<
        set -eux
        sniffles --input ~{sep=" " snfs} \
            --vcf multisample.vcf
        tree
    >>>

    output {
        File vcf = "~{prefix}.vcf"
    }

    Int cpus = 8
    Int disk_size = 3*ceil(size(snfs, "GB"))
                                                                                                                                                                                                                                                                                                                                                                                                                                   #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             46,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
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