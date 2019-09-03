version 1.0

# workflow to be updated (currently doesn't work too well)

workflow MT_old_flow {
    input {
        File corrected_bam
        File corrected_bai
        File remaining_bam
        File remaining_bai
        File ref_fasta
        String mt_chr_name
        String SM
        String ID
        String docker
    }

    call AssembleMT {
        input:
            corrected_bam = corrected_bam,
            corrected_bai = corrected_bai,
            remaining_bam = remaining_bam,
            remaining_bai = remaining_bai,
            ref_fasta = ref_fasta,
            mt_chr_name = mt_chr_name,
            SM = SM,
            ID = ID,
            docker = docker
    }
}

task AssembleMT {
    input {
        File corrected_bam
        File corrected_bai
        File remaining_bam
        File remaining_bai
        File ref_fasta
        String mt_chr_name
        String SM
        String ID
        String docker
    }

    Int cpus = 4
    Int disk_size = 4*ceil(size(corrected_bam, "GB") + size(corrected_bai, "GB") + size(remaining_bam, "GB") + size(remaining_bai, "GB") + size(ref_fasta, "GB"))

    String reads = "reads.fastq.gz"

    command <<<
        set -euxo pipefail

        # select MT
        samtools view -hb ~{corrected_bam} ~{mt_chr_name} | samtools fastq - | gzip -1 > ~{reads}
        samtools view -hb ~{remaining_bam} ~{mt_chr_name} | samtools fastq - | gzip -1 >> ~{reads}

        # align
        minimap2 -x ava-pb -t ~{cpus} ~{reads} ~{reads} | gzip -1 > reads.paf.gz

        # layout
        miniasm -f ~{reads} reads.paf.gz > reads.gfa
        awk '$1 ~/S/ { print "\x3E" $2 "\n" $3 }' reads.gfa > reads.fasta

        # correct 1
        minimap2 -t ~{cpus} reads.fasta ~{reads} > reads.gfa1.paf
        racon -t ~{cpus} ~{reads} reads.gfa1.paf reads.fasta > reads.racon1.fasta

        # correct 2
        minimap2 -t ~{cpus} reads.racon1.fasta ~{reads} > reads.gfa2.paf
        racon -t ~{cpus} ~{reads} reads.gfa2.paf reads.racon1.fasta > mt.fasta

        # align to ref
        minimap2 -ayY --MD --eqx -x asm20 -R '@RG\tID:~{ID}\tSM:~{SM}' -t ~{cpus} ~{ref_fasta} mt.fasta | samtools sort -@~{cpus} -m4G -o mt.bam
        samtools index mt.bam

        # call variants
        bcftools mpileup -Ou -f ~{ref_fasta} mt.bam | bcftools call -mv -Ov --ploidy 1 -o mt.vcf
    >>>

    output {
        File mt_asm = "mt.fasta"
        File mt_bam = "mt.bam"
        File mt_bai = "mt.bam.bai"
        File mt_vcf = "mt.vcf"
    }

    runtime {
        cpu: "~{cpus}"
        memory: "20G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}