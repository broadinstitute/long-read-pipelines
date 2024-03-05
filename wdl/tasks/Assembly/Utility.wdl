version 1.0

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

task mergefq{
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
        
    >>>

    output{
        File merged_fq = "merge.fq"
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

task hifiasm_asm{
    input{
        File reads
        String prefix
        Int num_cpus = 4
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
        memory: "64 GiB"
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

task ragtag_construct{
    input{
        File asm
        File ref
        String outputfolder
    }
    command <<<
    set -x pipefail
    ragtag.py scaffold ~{ref} ~{asm} -o ~{outputfolder}

    >>>

    output{
        File scaffold="~{outputfolder}/ragtag.scaffold.fasta"
        File stats = "~{outputfolder}/ragtag.scaffold.stats"
        File paf = "~{outputfolder}/ragtag.scaffold.asm.paf"
        File confidence = "~{outputfolder}/ragtag.scaffold.confidence.txt"
    }

    Int disk_size = 10 

    runtime {
        cpu: 2
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/ragtag:v1"
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
        #samtools index ~{pref}.~{region}.bam
    >>>

    output{
        File local_fa="~{pref}.fasta"
        #File local_bai="~{pref}.~{region}.bai"
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


task tsggapcloser_gapfilling{
    input{
        File scaffold_input
        File read_fq
        String outputprefix
    }
    command <<<
        set -x pipefail
        tgsgapcloser --scaff ~{scaffold_input} --reads ~{read_fq} --output ~{outputprefix} --ne --tgstype pb >pipe.log 2>pipe.err

    >>>

    output{
        File scaffold="~{outputprefix}.scaff_seqs"
        #File stats = "~{outputprefix}.gap_fill_detail"
    }

    Int disk_size = 10 

    runtime {
        cpu: 2
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/tgsgapcloser:v3"
    }
}