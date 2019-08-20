workflow TargetedAsm {
    Array[String] regions

    String parent1_corrected_bam
    String parent1_corrected_bai

    String parent2_corrected_bam
    String parent2_corrected_bai

    String child_corrected_bam
    String child_corrected_bai

    String ref_fasta
    String ref_fasta_fai
    String ref_dict

    String docker_image="kgarimella/pbasm@sha256:01661003c16177b6a2112a75f3ee9e4bd2ca41d0c77076b40ebf33f499a4356b"

    call SelectRegions as SelectParent1 {
        input:
            regions=regions,
            input_bam=parent1_corrected_bam,
            input_bai=parent1_corrected_bai,
            docker_image=docker_image
    }

    call SelectRegions as SelectParent2 {
        input:
            regions=regions,
            input_bam=parent2_corrected_bam,
            input_bai=parent2_corrected_bai,
            docker_image=docker_image
    }

    call SelectRegions as SelectChild {
        input:
            regions=regions,
            input_bam=child_corrected_bam,
            input_bai=child_corrected_bai,
            docker_image=docker_image
    }

    call BinRegions {
        input:
            reads_parent1=SelectParent1.reads,
            reads_parent2=SelectParent2.reads,
            reads_child=SelectChild.reads,
            docker_image="kgarimella/pbmetrics@sha256:4ce0a5d5172cbddca8f9512719f94b11749272b4ca8f1a54dd31713010ea37ab"
    }

    call Minimap2 as AlignReadsBin1 {
        input:
            reads=BinRegions.reads_p1,
            ref_fasta=ref_fasta,
            rgid="AlignReadsBin1",
            sample_name="child_readsbin1",
            docker_image=docker_image,
    }

    call Minimap2 as AlignReadsBin2 {
        input:
            reads=BinRegions.reads_p2,
            ref_fasta=ref_fasta,
            rgid="AlignReadsBin2",
            sample_name="child_readsbin2",
            docker_image=docker_image,
    }

    call AssembleRegions as AssembleBin1 {
        input:
            reads=BinRegions.reads_p1,
            docker_image=docker_image
    }

    call Minimap2 as AlignBin1 {
        input:
            reads=AssembleBin1.racon2,
            ref_fasta=ref_fasta,
            rgid="AssembleBin1",
            sample_name="child_bin1",
            docker_image=docker_image,
    }

    call AssembleRegions as AssembleBin2 {
        input:
            reads=BinRegions.reads_p2,
            docker_image=docker_image
    }

    call Minimap2 as AlignBin2 {
        input:
            reads=AssembleBin2.racon2,
            ref_fasta=ref_fasta,
            rgid="AssembleBin2",
            sample_name="child_bin2",
            docker_image=docker_image,
    }

    call AssembleRegions as AssembleParent1 {
        input:
            reads=SelectParent1.reads,
            docker_image=docker_image
    }

    call AssembleRegions as AssembleParent2 {
        input:
            reads=SelectParent2.reads,
            docker_image=docker_image
    }

    call AssembleRegions as AssembleChild {
        input:
            reads=SelectChild.reads,
            docker_image=docker_image
    }

    call Minimap2 as AlignParent1 {
        input:
            reads=AssembleParent1.racon2,
            ref_fasta=ref_fasta,
            rgid="AssembleParent1",
            sample_name="parent1",
            docker_image=docker_image,
    }

    call Minimap2 as AlignParent2 {
        input:
            reads=AssembleParent2.racon2,
            ref_fasta=ref_fasta,
            rgid="AssembleParent2",
            sample_name="parent2",
            docker_image=docker_image,
    }

    call Minimap2 as AlignChild {
        input:
            reads=AssembleChild.racon2,
            ref_fasta=ref_fasta,
            rgid="AssembleChild",
            sample_name="child",
            docker_image=docker_image,
    }
}

task SelectRegions {
    Array[String] regions
    File input_bam
    File input_bai
    String docker_image

    Int cpus = 2
    Int disk_size = ceil(2*(size(input_bam, "GB") + size(input_bai, "GB")))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        samtools view -hb ${input_bam} ${sep=' ' regions} | samtools fastq - | gzip -1 > reads.fastq.gz
        samtools view -f 0x4 ${input_bam} | samtools fastq - | gzip -1 >> reads.fastq.gz

        df -h .
        tree -h
    >>>

    output {
        File reads = "reads.fastq.gz"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "40G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
    }
}

task BinRegions {
    File reads_parent1
    File reads_parent2
    File reads_child
    String docker_image

    Int cpus = 2

    command <<<
        mccortex 21 build -f -m 20G -k 21 -S -s parent1 -1 ${reads_parent1} -s parent2 -1 ${reads_parent2} -s child -1 ${reads_child} trio_regions.ctx
        /bin_reads.py trio_regions.ctx ${reads_child} reads.parent1.fastq.gz reads.parent2.fastq.gz binning.stats.txt
    >>>

    output {
        File trio_regions = "trio_regions.ctx"
        File reads_p1 = "reads.parent1.fastq.gz"
        File reads_p2 = "reads.parent2.fastq.gz"
        File binning_stats = "binning.stats.txt"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "40G"
        bootDiskSizeGb: 20
        disks: "local-disk 40 SSD"
    }
}

task AssembleRegions {
    File reads
    String docker_image

    Int cpus = 8
    #Int disk_size = "10G"

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        # align
        minimap2 -x ava-pb -t ${cpus} ${reads} ${reads} | gzip -1 > reads.paf.gz

        # layout
        miniasm -f ${reads} reads.paf.gz > reads.gfa
        awk '$1 ~/S/ { print "\x3E" $2 "\n" $3 }' reads.gfa > reads.fasta

        # correct 1
        minimap2 -t ${cpus} reads.fasta ${reads} > reads.gfa1.paf
        racon -t ${cpus} ${reads} reads.gfa1.paf reads.fasta > reads.racon1.fasta

        # correct 2
        minimap2 -t ${cpus} reads.racon1.fasta ${reads} > reads.gfa2.paf
        racon -t ${cpus} ${reads} reads.gfa2.paf reads.racon1.fasta > reads.racon2.fasta

        df -h .
        tree -h
    >>>

    output {
        File paf = "reads.paf.gz"
        File gfa = "reads.gfa"
        File paf1 = "reads.gfa1.paf"
        File paf2 = "reads.gfa2.paf"
        File racon1 = "reads.racon1.fasta"
        File racon2 = "reads.racon2.fasta"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "40G"
        bootDiskSizeGb: 20
        disks: "local-disk 40 SSD"
    }
}

task Minimap2 {
    File ref_fasta
    File reads
    String rgid
    String sample_name
    String docker_image

    Int cpus = 3
    Int disk_size = ceil(size(ref_fasta, "GB")) + 4*ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        minimap2 -ayY --MD --eqx -x asm10 -R '@RG\tID:${rgid}\tSM:${sample_name}' -t ${cpus} ${ref_fasta} ${reads} | samtools sort -@${cpus} -m4G -o regions.aligned.sorted.bam
        samtools index regions.aligned.sorted.bam

        df -h .
        tree -h
    >>>

    output {
        File aligned_bam = "regions.aligned.sorted.bam"
        File aligned_bai = "regions.aligned.sorted.bam.bai"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        disks: "local-disk ${disk_size} SSD"
        memory: "20G"
        bootDiskSizeGb: 20
        preemptible: 3
        maxRetries: 3
    }
}
