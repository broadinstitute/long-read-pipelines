version 1.0

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
        String docker
    }

    Int cpus = 4
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

    runtime {
        cpu: "~{cpus}"
        memory: "20G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}