version 1.0

workflow Bam2Fastq{
    meta{
        description: "a workflow that extract vcfs and bams in a genomic interval"
    }
    input{
        File bam
        String prefix
    }
    call bam2fastq{input: bam_input=bam, pref=prefix}
    output{
        File fq1 = bam2fastq.fq1
        File fq2 = bam2fastq.fq2
    }
}



task bam2fastq{
    input{
        File bam_input
        String pref
    }
    command <<<
        samtools index ~{bam_input}
        bedtools bamtofastq -i ~{bam_input} \
                      -fq ~{pref}.aln.R1.fq \
                      -fq2 ~{pref}.aln.R2.fq
    >>>

    output{
        File fq1="~{pref}.aln.R1.fq"
        File fq2 = "~{pref}.aln.R2.fq"
    }

    Int disk_size = 10 + ceil(2 * size(bam_input, "GiB"))

    runtime {
        cpu: 1
        memory: "6 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}
