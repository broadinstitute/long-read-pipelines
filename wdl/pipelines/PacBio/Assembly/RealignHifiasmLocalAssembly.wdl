version 1.0

workflow RealignHifiasmLocalAssembly{
    input{
        File wholegenomebam
        File wholegenomebai
        String asmregion
        String prefix
        Int nthreads
    }
    call extract_reads{input: bam_input=wholegenomebam, bam_index=wholegenomebai, region=asmregion, pref=prefix}
    call hifiasm_asm{input: reads=extract_reads.local_fq, prefix=prefix}
    call realign_reads{input: wholebam=wholegenomebam, wholebai=wholegenomebai ,assembly=hifiasm_asm.assembly_primary, num_threads=nthreads, pref=prefix}
    call hifiasm_asm_reassemble{input: first_round_reads=extract_reads.local_fq, reads=realign_reads.realignedreadsfq, prefix=prefix}

    meta{
        Purpose:"Local assembly using hifiasm"
    }
    output{
        File realignedfq=realign_reads.realignedreadsfq
        File first_assembly_hap1=hifiasm_asm.assembly_hap1
        File first_assembly_hap2=hifiasm_asm.assembly_hap2
        File second_assembly_hap1=hifiasm_asm_reassemble.assembly_hap1
        File second_assembly_hap2=hifiasm_asm_reassemble.assembly_hap2
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
        #samtools fastq ~{pref}.bam > ~{pref}.fastq
        bedtools bamtofastq  -i ~{pref}.bam -fq ~{pref}.fastq
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
        File assembly_primary="~{prefix}.bp.p_ctg.fa"
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

task hifiasm_asm_reassemble{
    input{
        File first_round_reads
        File reads
        String prefix
        Int num_cpus = 32
    }

    Int disk_size = 1 + ceil(2 * size(reads, "GiB"))

    command <<<

        set -euxo pipefail
        cat ~{first_round_reads} ~{reads} > merge.fq
        hifiasm -o ~{prefix} -t ~{num_cpus} merge.fq
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