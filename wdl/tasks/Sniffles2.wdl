version 1.0


import "Structs.wdl"

# "This workflow calls SV candidates using Sniffles2 population mode.
# 1) Call SV candidates and create .snf file for each sample
#  2) Combined calling using multiple .snf files into a single .vcf"

workflow Sniffles2 {

    input {
        Array[File] sampleBAMs
        Array[File] sampleBai
        Int minsvlen
        String prefix
        String sample_id}

    scatter (pair in zip(sampleBAMs, sampleBai)) {
        call SampleSV {
                input:
                    bam = pair.left,
                    bai = pair.right,
                    minsvlen = minsvlen,
                    prefix = prefix,
                    sample_id = sample_id
                }
         }

    call MergeCall { input: snfs = SampleSV.snf}


    output {
         Array[File] single_snf = SampleSV.snf
         File multisample_vcf = MergeCall.vcf
        }
}



task SampleSV {

    input {
        File bam
        File bai
        Int minsvlen
        String? chr
        String sample_id
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:              "input BAM from which to call SVs"
        bai:              "index accompanying the BAM"
        minsvlen:         "minimum SV length in bp. Default 35"
        chr:              "chr on which to call variants"
        sample_id:        "Sample ID"
        prefix:           "prefix for output"
    }

    Int cpus = 8
    Int disk_size = 2*ceil(size([bam, bai], "GB"))
    String fileoutput = if defined(chr) then "~{prefix}.~{chr}.sniffles.snf" else "~{prefix}.sniffles.snf"
    String vcf_output = if defined(chr) then "~{prefix}.~{chr}.sniffles.vcf" else "~{prefix}.sniffles.vcf"

    command <<<
        set -x

        sniffles -t ~{cpus} \
                 -i ~{bam} \
                 --minsvlen ~{minsvlen} \
                 --sample-id ~{sample_id} \
                 --vcf ~{vcf_output} \
                 --snf ~{fileoutput}
        tree
    >>>

    output {
        File snf = "~{fileoutput}"
        File vcf = "~{vcf_output}"
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


task MergeCall {
    input {
        Array[File] snfs
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        snfs: ".snf files"
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