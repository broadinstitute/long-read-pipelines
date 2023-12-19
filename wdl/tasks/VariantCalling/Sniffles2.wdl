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
        String? output_bucket
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:              "input BAM from which to call SVs"
        bai:              "index accompanying the BAM"
        minsvlen:         "minimum SV length in bp. Default 50"
        sample_id:        "Sample ID"
        prefix:           "prefix for output"
        output_bucket:    "cloud path for output files"
    }

    Int disk_size = 2*ceil(size([bam, bai], "GB"))

    command <<<
        set -eux

        vcfName="~{prefix}.sniffles.vcf.gz"
        tbiName="~{prefix}.sniffles.vcf.gz.tbi"
        snfName="~{prefix}.sniffles.snf"
        sniffles -t "$(nproc)" \
                 -i "~{bam}" \
                 --minsvlen "~{minsvlen}" \
                 --sample-id "~{sample_id}" \
                 --vcf "${vcfName}" \
                 --snf "${snfName}"

        if ~{defined(output_bucket)}; then
            outDir=$(echo "~{output_bucket}" | sed 's+/?$+/+')
            gcloud storage cp "$vcfName" "$tbiName" "$snfName" "$outDir"
        fi
    >>>

    output {
        File snf = "$snfName"
        File vcf = "$vcfName"
        File tbi = "$tbiName"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             46,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sniffles2:2.2.1"
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
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sniffles2:2.2.1"
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
