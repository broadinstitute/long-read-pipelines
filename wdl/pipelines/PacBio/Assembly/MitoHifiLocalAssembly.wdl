version 1.0

workflow MitoHifiLocalAssembly{
    input{
        File ref_fa
        File ref_gb
        File baminput
        String region
        Int n_cpus
        String prefix
    }
    call extract_reads{input: bam_input=baminput, region=region, pref=prefix}
    call MitoHifiAsm{input: reads=extract_reads.local_fq, reffa=ref_fa, refgb=ref_gb, num_cpus=n_cpus}
    meta{
        Purpose:"Mitochondrial genome Local assembly using MitoHifi"
    }
    output{
        File MT_fa=MitoHifiAsm.final_fa
        File MT_gb=MitoHifiAsm.final_gb
        File MT_annotation_fig=MitoHifiAsm.final_annotation_fig
        File MT_coverage_fig=MitoHifiAsm.final_coverage_fig
        File MT_stats=MitoHifiAsm.final_stats 
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

task MitoHifiAsm{
    input{
        File reads
        File reffa
        File refgb
        Int num_cpus
    }

    Int disk_size = 50 + ceil(2 * size(reads, "GiB"))

    command <<<

        set -euxo pipefail

        mitohifi.py -r ~{reads} -f ~{reffa} -g ~{refgb} -t ~{num_cpus} -o 1 
    >>>

    output{
        File final_fa="final_mitogenome.fasta"
        File final_gb="final_mitogenome.gb"
        File final_annotation_fig="final_mitogenome.annotation.png"
        File final_coverage_fig="final_mitogenome.coverage.png"
        File final_stats="contigs_stats.tsv"

    }
    runtime {
        cpu: num_cpus
        memory: "10 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        #docker: "ghcr.io/marcelauliano/mitohifi:master"
        docker:"biocontainers/mitohifi:3.0.0_cv1"
    }    
}