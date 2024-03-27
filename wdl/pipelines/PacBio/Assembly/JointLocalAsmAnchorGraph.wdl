version 1.0

import "../../../tasks/Utility/Utils.wdl" as U

workflow JointLocalAsmAnchorGraph{
    meta{
        description: "a workflow that extract vcfs and bams in a genomic interval"
    }
    input{
        Array[File] wholegenomebam
        Array[File] wholegenomebai
        String genomeregion
        String prefix
    }
    scatter (bam_bai in zip(wholegenomebam, wholegenomebai)){
        File bam = bam_bai.left
        File bai = bam_bai.right
        call U.InferSampleName { input: 
                    bam = bam, 
                    bai = bai
                }
        String sample_id = InferSampleName.sample_name
        call U.SubsetBam as Sub {input:
            bam = bam,
            bai = bai,
            locus = genomeregion,
            prefix = sample_id
        }
        call extract_bam_addsample{input: bam_input=Sub.subset_bam, bam_index=Sub.subset_bai, region= genomeregion, pref=sample_id}
    }

    call catfasta{input: fas= extract_bam_addsample.local_fa, pref=prefix}
    
    output{
        File merged_fa = catfasta.fasta
        #File subset_bai = extract_bam.local_bai
    }
}



task extract_bam_addsample{
    input{
        File bam_input
        File bam_index
        String region
        String pref
    }
    command <<<
        #samtools view --with-header ~{bam_input} -b ~{region} -o ~{pref}.bam
        samtools fasta ~{bam_input} > ~{pref}.fasta
        #bedtools bamtofastq -i ~{pref}.bam -fq ~{pref}.fastq

        #sed -n '1~4s/^@/>/p;2~4p' ~{pref}.fastq > ~{pref}.fasta
        sed '/^>/s/$/\|~{pref}/' < ~{pref}.fasta > ~{pref}_out.fasta # add sample information to each of the fasta file
        
    >>>

    output{
        File local_fa="~{pref}_out.fasta"
        # File local_bam = "~{pref}.bam"
        # File local_bai="~{pref}.~{region}.bai"
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

task catfasta{
    input{
        Array[File] fas
        String pref
    }
    command <<<
        cat  ~{sep=" " fas} > ~{pref}.aln.fa
    >>>

    output{
        File fasta="~{pref}.aln.fa"
    }

    Int disk_size = 100

    runtime {
        cpu: 1
        memory: "6 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}