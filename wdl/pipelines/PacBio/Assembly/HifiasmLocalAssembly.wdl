version 1.0

workflow HifiasmLocalAssembly{
    input{
        File baminput
        String region
        Int n_cpus
        String prefix
    }
    call extract_reads{input: bam_input=baminput, region=region, pref=prefix}
    call hifiasm_asm{input: reads=extract_reads.local_fq, prefix=prefix}
    meta{
        Purpose:"Local assembly using hifiasm"
    }
}
task extract_reads{
    input{
        File bam_input
        String region
        String pref
    }
    command <<<
        samtools index ~{bam_input}
        samtools view --with-header ~{bam_input} -b ~{region} -o ~{pref}.bam
        samtools fastq ~{pref}.bam > ~{pref}.fastq
    >>>

    output{
        File local_fq="~{pref}.fastq"
    }

    Int disk_size = 1 + ceil(2 * size(bam_input, "GiB"))

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

task hifiasm_asm{
    input{
        File reads
        String prefix
        Int num_cpus = 32
    }

    Int disk_size = 1 + ceil(2 * size(reads, "GiB"))

    command <<<

        set -euxo pipefail

        hifiasm -o ~{prefix} -t~{num_cpus} ~{reads}
        awk '/^S/{print ">"$2; print $3}' ~{prefix}.bp.p_ctg.gfa > ~{prefix}.bp.p_ctg.fa
        awk '/^S/{print ">"$2;print $3}' ~{prefix}.bp.hap1.p_ctg.gfa > ~{prefix}.bp.hap1.p_ctg.fa
        awk '/^S/{print ">"$2;print $3}' ~{prefix}.bp.hap2.p_ctg.gfa > ~{prefix}.bp.hap2.p_ctg.fa
    >>>

    output{
        File assembly_hap1="~{prefix}.bp.hap1.p_ctg.fa"
        File assembly_hap2="~{prefix}.bp.hap2.p_ctg.fa"
    }
    runtime {
        cpu: num_cpus
        memory: "10 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/assembly:v1"
    }    
}