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

    parameter_meta {
        bam:              { desciption: "input BAM from which to call SVs", localization_optional: true }
        bai:              "index accompanying the BAM"
        minsvlen:         "minimum SV length in bp. Default 50"
        sample_id:        "Sample ID"
        prefix:           "prefix for output"
        phase_sv:         "if you're sure the BAM is phased/haplotagged, turn this on to generate phased SV"
        tandem_repeat_bed: "BED file containing TRF finder (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"
    }

    input {
        File bam
        File bai
        Int minsvlen = 50
        String sample_id
        String prefix
        File? tandem_repeat_bed
        Boolean phase_sv = false
        RuntimeAttr? runtime_attr_override
    }

    String postfix = if phase_sv then "-phased" else ""
    String snf_output = "~{prefix}.sniffles~{postfix}.snf"
    String vcf_output = "~{prefix}.sniffles~{postfix}.vcf.gz"
    String tbi_output = "~{prefix}.sniffles~{postfix}.vcf.gz.tbi"

    String local_bam = "/cromwell_root/~{basename(bam)}"

    command <<<
        set -euxo pipefail

        time gcloud storage cp ~{bam} ~{local_bam}
        mv ~{bai} "~{local_bam}.bai"

        touch ~{bai}  # handle the bai-older-than-bam warning
        sniffles -t ~{cpus} \
                 -i ~{local_bam} \
                 --minsvlen ~{minsvlen} \
                 --sample-id ~{sample_id} \
                 ~{if defined(tandem_repeat_bed) then "--tandem-repeats ~{tandem_repeat_bed}" else ""} \
                 ~{true="--phase" false="" phase_sv} \
                 --vcf ~{vcf_output} \
                 --snf ~{snf_output}
    >>>

    output {
        File snf = "~{snf_output}"
        File vcf = "~{vcf_output}"
        File tbi = "~{tbi_output}"
    }

    #########################

    Int cpus = 8
    Int disk_size = 2*ceil(size([bam, bai], "GB"))

    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sniffles2:2.2"
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

    #########################
    Int cpus = 8
    Int disk_size = 3*ceil(size(snfs, "GB"))

    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             46,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sniffles2:2.2"
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
