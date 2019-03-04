
workflow RealignBamWorkflow {
    File ref_fasta
    File input_bam
    File input_bam_bai

    call RealignBam {
        input:
            ref_fasta=ref_fasta,
            input_bam=input_bam,
            input_bam_bai=input_bam_bai
    }
}

task RealignBam {
    File ref_fasta
    File input_bam
    File input_bam_bai
    Int disk_size = ceil(size(ref_fasta, "GB") + 2.5*size(input_bam, "GB") + 10)

    File ref_fasta_fai = "${ref_fasta}.fai"
    File ref_fasta_amb = "${ref_fasta}.amb"
    File ref_fasta_ann = "${ref_fasta}.ann"
    File ref_fasta_bwt = "${ref_fasta}.bwt"
    File ref_fasta_genome = "${ref_fasta}.genome"
    File ref_fasta_pac = "${ref_fasta}.pac"
    File ref_fasta_sa = "${ref_fasta}.sa"


    command {
        set -xe

        ls -lh ${input_bam}
        wget -q https://github.com/ssadedin/bazam/releases/download/1.0.1/bazam.jar
        wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
        tar xjf bwa-0.7.17.tar.bz2
        cd bwa-0.7.17
        make
        mv bwa ../
        cd ..
        java -Xmx20G -jar ./bazam.jar -bam ${input_bam} | ./bwa mem -p ${ref_fasta} - | samtools view -bSu -  | samtools sort -o out.bam
        samtools index out.bam
    }

    output {
        File fasta_index="${ref_fasta}*"
        File output_bam="out.bam"
        File output_bai="out.*bai"
    }

    runtime {
        docker: "weisburd/base-image-for-str-tools@sha256:bc02c67c69bbd13165ef2cd83f7b2ec87814d99fc007070cc7cd8865b29998a9"
        disks: "local-disk ${disk_size} HDD"
        memory: "32 GB"
        cpu: "4"
    }
}
