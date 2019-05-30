workflow MHCAsmWorkflow {
    String input_bam
    String input_bai
    String docker_image="kgarimella/pbasm@sha256:01661003c16177b6a2112a75f3ee9e4bd2ca41d0c77076b40ebf33f499a4356b"

    call AssembleMHC {
        input:
            input_bam=input_bam,
            input_bai=input_bai,
            docker_image=docker_image
    }
}

task AssembleMHC {
    File input_bam
    File input_bai
    String docker_image

    Int cpus = 8
    Int disk_size = ceil(2*(size(input_bam, "GB")))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        #export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

        # extract
        samtools view -h ${input_bam} chr6:28000000-34000000 | samtools fastq - > reads.fq

        # align
        minimap2 -x ava-pb -t ${cpus} reads.fq reads.fq | gzip -1 > reads.paf.gz

        # layout
        miniasm -f reads.fq reads.paf.gz > reads.gfa
        awk '$1 ~/S/ { print ">" $2 "\n" $3 }' reads.gfa > reads.fasta

        # correct 1
        minimap2 -t ${cpus} reads.fasta reads.fq > reads.gfa1.paf
        racon -t ${cpus} reads.fq reads.gfa1.paf reads.fasta > reads.racon1.fasta

        # correct 2
        minimap2 -t ${cpus} reads.racon1.fasta reads.fq > reads.gfa2.paf
        racon -t ${cpus} reads.fq reads.gfa2.paf reads.racon1.fasta > reads.racon2.fasta

        df -h .
        tree -h
    >>>

    output {
        File racon1 = "reads.racon1.fasta"
        File racon2 = "reads.racon2.fasta"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "40G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
    }
}
