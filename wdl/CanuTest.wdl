version 1.0

workflow CanuTest {
    input {
        File bam
        File bai
        String platform
        Boolean is_corrected
    }

    String docker_asm = "kgarimella/lr-asm:0.01.09"

    call CanuTarget {
        input:
            bam=bam,
            bai=bai,
            region="chrM",
            platform=platform,
            is_corrected=is_corrected,
            genome_size="16.6k",
            prefix="mt",
            docker=docker_asm
    }
}

task CanuTarget {
    input {
        File bam
        File bai
        String region
        String platform
        Boolean is_corrected
        String genome_size
        String prefix
        String docker
    }

    String data_type = "-" + (if platform == "PACBIO" then "pacbio" else "nanopore") + "-" + (if is_corrected then "corrected" else "raw")
    Int disk_size = 4*ceil(size(bam, "GB") + size(bai, "GB"))

    command <<<
        set -euxo pipefail

        # select region
        samtools view -hb ~{bam} ~{region} | samtools fastq - | gzip -1 > reads.fastq.gz

        # assemble MT
        canu -p ~{prefix} -d ./ genomeSize=~{genome_size} readSamplingCoverage=500 ~{data_type} reads.fastq.gz
    >>>

    output {
        File report          = "mt.report"

        File contigs_fasta   = "mt.contigs.fasta"
        File unassembled     = "mt.unassembled.fasta"
        File unitigs_fasta   = "mt.unitigs.fasta"

        File contigs_layout  = "mt.contigs.layout"
        File unitigs_layout  = "mt.unitigs.layout"
        File unitigs_bed     = "mt.unitigs.bed"

        File contigs_gfa     = "mt.contigs.gfa"
        File unitigs_gfa     = "mt.unitigs.gfa"

        File corrected_reads = "mt.correctedReads.fasta.gz"
        File trimmed_reads   = "mt.trimmedReads.fasta.gz"
    }

    runtime {
        cpu: 8
        memory: "12G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}
