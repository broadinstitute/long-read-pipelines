version 1.0

workflow AlignAsmtoGenome{
    meta{
        description: "a workflow that extract multiple sample in a genomic interval and construct a graph"
    }
    input{
        File assembly
        File reference
        String prefix
    }

    
    call minimap_align{input: query_fasta=assembly, reference_fasta= reference, pref=prefix}
    


    
    output{
        # File merged_fa = catfasta.fasta
        File bam_file = minimap_align.bam
        File bai_file = minimap_align.bai
    }
}



task minimap_align{
    input{
        File query_fasta
        File reference_fasta
        String pref
    }
    command <<<
    minimap2 -ax asm5 ~{reference_fasta} ~{query_fasta} > ~{pref}.aln.sam

    samtools view -H ~{pref}.aln.sam

    samtools sort ~{pref}.aln.sam -o ~{pref}.aln.sorted.bam
    samtools index ~{pref}.aln.sorted.bam ~{pref}.aln.sorted.bam.bai

    >>>

    output{
        File bam="~{pref}.aln.sorted.bam"
        File bai= "~{pref}.aln.sorted.bam.bai"
    }

    Int disk_size = 100

    runtime {
        cpu: 1
        memory: "128 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
    }
}