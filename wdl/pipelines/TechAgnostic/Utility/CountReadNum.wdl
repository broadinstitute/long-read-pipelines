version 1.0

workflow ExtractVcfandBam{
    meta{
        description: "a workflow that extract vcfs and bams in a genomic interval"
    }
    input{
        File bam
        File bai
        String genomeregion
        String prefix
    }
    call count_read_num{input: bam_input=bam, bam_index=bai, region= genomeregion, pref=prefix}
    output{
        Int read_number = read_int(count_read_num.read_number) 
    }
}

task count_read_num{
    input{
        File bam_input
        File bam_index
        String region
        String pref
    }
    command <<<
        samtools view -c ~{bam_input} > read_num.txt
        #samtools index ~{pref}.~{region}.bam
    >>>

    output{
        File read_number = "read_num.txt"
    }

    Int disk_size = 10 + ceil(2 * size(bam_input, "GiB"))

    runtime {
        cpu: 1
        memory: "10 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}