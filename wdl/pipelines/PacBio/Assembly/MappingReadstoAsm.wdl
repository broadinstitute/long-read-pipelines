version 1.0

workflow MappingReadstoAsm{
    input{
        File baminput
        File baiinput
        File asm
        String outputpref
        Int num_cpus
    }
    call realign_reads{input: wholebam=baminput, wholebai=baiinput, assembly=asm, num_threads = num_cpus, pref=outputpref}
    #call whatshap_phasing{input: inputbam=baminput, ref=reference, outputprefix=outputpref, vcf=bcf_mpileup.bcf_file}
    meta{
        Purpose:"Mapping reads to asm to attract more reads"
    }
    output{
        
        File realignedbam=realign_reads.realignedbam

    }
}
task realign_reads{
    input{
        File wholebam
        File wholebai
        File assembly
        String num_threads
        String pref
    }
    command <<<
        set -euxo pipefail
        # bedtools bamtofastq  -i ~{wholebam} -fq ~{pref}.whole.fastq
        samtools fastq -f 4 ~{wholebam} | minimap2 -t ~{num_threads} -ax map-hifi ~{assembly} - | samtools sort -m4G -@4 -o ~{pref}_realign.sorted.bam - 
        samtools index ~{pref}_realign.sorted.bam
        samtools fastq -F 4 ~{pref}_realign.sorted.bam > ~{pref}_realign.fastq # extract aligned reads
        #rm ~{pref}.sam
    >>>

    output{
        File realignedreadsfq= "~{pref}_realign.fastq"
        File realignedbam="~{pref}_realign.sorted.bam"
        #File alignmentbam="~{pref}.sorted.bam"
        #File alignmentbai="~{pref}.sorted.bam.bai"
    }
    Int disk_size = 100 + ceil(2 * (size(wholebam, "GiB") + size(assembly, "GiB") ))

    runtime {
        cpu: 1
        memory: "20 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
    }
}