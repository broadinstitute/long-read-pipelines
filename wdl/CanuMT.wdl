version 1.0

import "Structs.wdl"

# TODO: describe purpose
task CanuMT {
    input {
        File corrected_bam
        File corrected_bai
        File remaining_bam
        File remaining_bai
        File ref_fasta
        String mt_chr_name
        String SM
        String ID

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(corrected_bam, "GB") + size(corrected_bai, "GB") + size(remaining_bam, "GB") + size(remaining_bai, "GB") + size(ref_fasta, "GB"))

    String creads = "creads.fastq.gz"
    String rreads = "rreads.fastq.gz"

    command <<<
        set -euxo pipefail

        # select MT
        samtools view -hb ~{corrected_bam} ~{mt_chr_name} | samtools fastq - | gzip -1 > ~{creads}
        samtools view -hb ~{remaining_bam} ~{mt_chr_name} | samtools fastq - | gzip -1 > ~{rreads}

        canu -p mt -d mttest genomeSize=16.6k -pacbio-raw ~{rreads}

        find . -type f -exec ls -lh {} \;
    >>>

    output {
        Array[File] all = glob("mttest/*")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4, 
        mem_gb:             20, 
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "kgarimella/lr-asm:0.01.06"
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