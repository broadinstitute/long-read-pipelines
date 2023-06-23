version 1.0
workflow SvimAsmAssemblyVariants{
    input{
        File ref_fa
        File assembly_fa
        Int n_cpus
        String prefix
    }
    call alignreads{input: assembly=assembly_fa, reference=ref_fa, num_threads=n_cpus, pref=prefix}
    call CallVariants{input: bamfile=alignreads.alignmentbam, baifile=alignreads.alignmentbai, reference=ref_fa}
    meta{
        Purpose:"Call Variants from assembly using SvimAsm"
    }
    output{
        Array[File] Final=CallVariants.a_files
        File vcf=CallVariants.vcf
        # File MT_gb=MitoHifiAsm.final_gb
        # File MT_annotation_fig=MitoHifiAsm.final_annotation_fig
        # File MT_coverage_fig=MitoHifiAsm.final_coverage_fig
        # File MT_stats=MitoHifiAsm.final_stats 
    }
}

task alignreads{
    input{
        File assembly
        File reference
        String num_threads
        String pref
    }
    command <<<
        minimap2 -a -x asm5 --cs -r2k -t ~{num_threads} ~{reference} ~{assembly} > ~{pref}.sam
        samtools sort -m4G -@4 -o ~{pref}.sorted.bam ~{pref}.sam
        samtools index ~{pref}.sorted.bam
    >>>

    output{
        File alignmentbam="~{pref}.sorted.bam"
        File alignmentbai="~{pref}.sorted.bam.bai"
    }
    Int disk_size = 1 + ceil(2 * (size(assembly, "GiB") + size(reference, "GiB") ))

    runtime {
        cpu: 1
        memory: "10 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
    }
}

task CallVariants{
    input{
        File bamfile
        File baifile
        File reference

    }
    command <<<
        svim-asm haploid . ~{bamfile} ~{reference}
    >>>

    output{
        Array[File] a_files = glob("*")
        File vcf="variants.vcf"
    }
    Int disk_size = 1 + ceil(2 * (size(bamfile, "GiB") ))

    runtime {
        cpu: 1
        memory: "10 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/svim-asm:v1"
    }
}