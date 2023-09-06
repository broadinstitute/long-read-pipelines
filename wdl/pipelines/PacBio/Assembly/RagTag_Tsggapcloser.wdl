version 1.0

workflow RagTag_Tsggapcloser{
    meta{
        description: "a workflow that use ragtag and tsggapcloser to scaffold a genome"
    }
    input{
        File input_fastas_hap1
        File input_fastas_hap2
        File reference_fasta
        File wholegenomebam
        File wholegenomebai
        String genomeregion
        String outputprefix
    }
    call ragtag_construct as ragtag_hap1{input: asm = input_fastas_hap1, ref = reference_fasta, outputfolder = outputprefix}
    call ragtag_construct as ragtag_hap2{input: asm = input_fastas_hap2, ref = reference_fasta, outputfolder = outputprefix}
    call extract_bam{input: bam_input=wholegenomebam, bam_index=wholegenomebai, region= genomeregion, pref=outputprefix}
    call tsggapcloser_gapfilling as tgs_hap1{input: scaffold_input = ragtag_hap1.scaffold, read_fq = extract_bam.local_fa, outputprefix = outputprefix}
    call tsggapcloser_gapfilling as tgs_hap2{input: scaffold_input = ragtag_hap2.scaffold, read_fq = extract_bam.local_fa, outputprefix = outputprefix}
    
    output{
        File ragtag_scaffold_hap1 = ragtag_hap1.scaffold
        File ragtag_scaffold_hap2 = ragtag_hap2.scaffold
        File tgs_scaffold_hap1 = tgs_hap1.scaffold
        File tgs_scaffold_hap2 = tgs_hap2.scaffold

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