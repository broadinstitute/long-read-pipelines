version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.14/wdl/tasks/Structs.wdl"
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.14/wdl/tasks/Utils.wdl" as Utils

workflow CallSmallVariants {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai
        File ref_dict
    }

    call Utils.MakeChrIntervalList { input: ref_dict = ref_dict }

    scatter (chr_info in MakeChrIntervalList.chrs) {
        call Longshot {
            input:
                bam = bam,
                bai = bai,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                chr = chr_info[0]
        }
    }

    call MergeLongshotCalls {
        input:
            vcfs = Longshot.vcf,
            ref_dict = ref_dict,
            prefix = basename(bam, ".bam")
    }

    output {
        File longshot_vcf = MergeLongshotCalls.vcf
        File longshot_tbi = MergeLongshotCalls.tbi
    }
}

task Longshot {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        String chr

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 4
    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB") + size(ref_fasta, "GB") + size(ref_fasta_fai, "GB"))
    String prefix = basename(bam, ".bam")

    command <<<
        set -euxo pipefail

        touch ~{prefix}.longshot.~{chr}.vcf
        longshot -F -r ~{chr} --bam ~{bam} --ref ~{ref_fasta} --out ~{prefix}.longshot.~{chr}.vcf
    >>>

    output {
        File vcf = "~{prefix}.longshot.~{chr}.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longshot:0.1.1"
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

task MergeLongshotCalls {
    input {
        Array[File] vcfs
        File ref_dict
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(vcfs, "GB")) + 1

    command <<<
        set -x

        VCF_WITH_HEADER=`ls -S ~{sep=' ' vcfs} | head -1`

        grep '^#' $VCF_WITH_HEADER | grep -v CHROM > header
        grep '^@SQ' ~{ref_dict} | awk '{ print "##contig=<ID=" $2 ",length=" $3 ">" }' | sed 's/[SL]N://g' >> header
        grep -m1 CHROM $VCF_WITH_HEADER >> header

        ((cat header) && (grep -h -v '^#' ~{sep=' ' vcfs})) | bcftools sort | bgzip > ~{prefix}.longshot.vcf.gz
        tabix -p vcf ~{prefix}.longshot.vcf.gz
    >>>

    output {
        File vcf = "~{prefix}.longshot.vcf.gz"
        File tbi = "~{prefix}.longshot.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             24,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longshot:0.1.1"
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
