version 1.0

#import "tasks/Utils.wdl" as Utils
#import "../../../tasks/Assembly/Canu.wdl" as Canu

workflow HiCanuLocalAssembly{
    meta{
        description: "a workflow that assemble the contigs using low coverage data"
    }
    input{
        File wholegenomebam
        File wholegenomebai
        String asmregion
        String prefix
        Int nthreads
        String region_size
    }
    call extract_reads{input: bam_input=wholegenomebam, bam_index=wholegenomebai, region=asmregion, pref=prefix}
    call hicanu_asm{input: reads=extract_reads.local_fq, prefix=prefix, genome_size=region_size}
    output{
        File canu_asm = hicanu_asm.canu_contigs_fasta
        
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

task hicanu_asm{
    input{
        File reads
        String prefix
        String genome_size
        
        Int num_cpus = 32
    }

    Int disk_size = 150 * ceil(size(reads, "GB"))

    command <<<

        set -euxo pipefail
        canu -p ~{prefix} -d . genomeSize=~{genome_size} minInputCoverage=3 stopOnLowCoverage=3 -pacbio-hifi ~{reads}
        #canu -p ~{prefix} -d . genomeSize=~{genome_size} minInputCoverage=3 stopOnLowCoverage=3 -pacbio-hifi ~{reads}
            >>>

    output{
        File canu_contigs_fasta = "~{prefix}.contigs.fasta"
    }
    runtime {
        cpu: num_cpus
        memory: "10 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-canu:0.1.0"
    }    
}

