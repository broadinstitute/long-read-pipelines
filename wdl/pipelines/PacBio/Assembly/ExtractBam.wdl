version 1.0

workflow ExtractBam{
    meta{
        description: "a workflow that extract vcfs and bams in a genomic interval"
    }
    input{
        File wholegenomebam
        File wholegenomebai
        String genomeregion
        String prefix
    }
    call extract_bam{input: bam_input=wholegenomebam, bam_index=wholegenomebai, region= genomeregion, pref=prefix}
    output{
        File subset_fa = extract_bam.local_fa
        File subset_bam = extract_bam.local_bam
        #File subset_bai = extract_bam.local_bai
    }
}



task extract_bam{
    input{
        File bam_input
        File bam_index
        String region
        String pref
    }
    command <<<
        samtools view --with-header ~{bam_input} -b ~{region} -o ~{pref}.bam
        samtools fasta ~{pref}.bam > ~{pref}.fasta
        sed '/^>/s/$/\_~{pref}/' < ~{pref}.fasta > ~{pref}_out.fasta # add sample information to each of the fasta file
        #samtools index ~{pref}.~{region}.bam
    >>>

    output{
        File local_fa="~{pref}_out.fasta"
        File local_bam = "~{pref}.bam"
        # File local_bai="~{pref}.~{region}.bai"
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
