version 1.0

workflow FlyeLocalAssembly{
    meta{
        description: "a workflow that assemble the contigs using low coverage data"
    }
    input{
        File wholegenomebam
        File wholegenomebai
        String asmregion
        String prefix
        Int nthreads
    }
    call extract_reads{input: bam_input=wholegenomebam, bam_index=wholegenomebai, region=asmregion, pref=prefix}
    call flye_asm{input: reads=extract_reads.local_fq, prefix=prefix}
    output{
        File flye_contig = flye_asm.contigs_fasta
        
    }
}

task extract_reads{
    input{
        File bam_input
        File bam_index
        String region
        String pref
    }
    command <<<
        samtools index ~{bam_input}
        samtools view --with-header ~{bam_input} -b ~{region} -o ~{pref}.bam
        samtools fastq ~{pref}.bam > ~{pref}.fastq
        #bedtools bamtofastq  -i ~{pref}.bam -fq ~{pref}.fastq
    >>>

    output{
        File local_fq="~{pref}.fastq"
    }

    Int disk_size = 100 + ceil(2 * size(bam_input, "GiB"))

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

task flye_asm{
    input{
        File reads
        String prefix
        Int num_cpus = 32
    }

    Int disk_size = 150 * ceil(size(reads, "GB"))

    command <<<

        set -euxo pipefail
        flye --pacbio-hifi ~{reads} --out-dir . --threads ~{num_cpus}
        >>>

    output{
        File contigs_fasta = "assembly.fasta"
    }
    runtime {
        cpu: num_cpus
        memory: "10 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-flye:2.8.3"
    }    
}

