version 1.0

workflow Bam2Fastq{
    meta{
        description: "a workflow that extract vcfs and bams in a genomic interval"
    }
    input{
        Array[File] fastas
        String prefix
    }
    call catfasta{input: fas=fastas, pref=prefix}
    output{
        File fq1 = catfasta.fasta
    }
}



task catfasta{
    input{
        Array[File] fas
        String pref
    }
    command <<<
        cat  ~{sep=" " fas} > ~{pref}.aln.fa
    >>>

    output{
        File fasta="~{pref}.aln.fa"
    }

    Int disk_size = 100

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
