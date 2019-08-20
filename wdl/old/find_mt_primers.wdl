workflow FindMTPrimersWorkflow {
    String input_bam
    String sm
    String docker_image="kgarimella/pbtools@sha256:f2775bd706bb1d955a2a77b915aed659b3ad68a18c61dd792d7729c2f4fdd4d9"

    call FindMTPrimers {
        input:
            input_bam=input_bam,
            sample_name=sm,
            docker_image=docker_image
    }
}

task FindMTPrimers {
    File input_bam
    String sample_name
    String docker_image

    Int cpus = 1
    Int disk_size = ceil(2*size(input_bam, "GB"))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        echo -e '>MethodA_primer1_2120F\nGGACACTAGGAAAAAACCTTGTAGAGAGAG\n>MethodA_primer2_2119R\nAAAGAGCTGTTCCTCTTTGGACTAACA\n>MethodB_primer1_16426F\nCCGCACAAGAGTGCTACTCTCCTC\n>MethodB_primer2_16425R\nGATATTGATTTCACGGAGGATGGTG' > mt_primers.fa
        samtools fastq ${input_bam} > input.fastq
        minimap2 -a input.fastq mt_primers.fa > ${sample_name}.mt_primers.sam

        df -h .
        tree -h
    >>>

    output {
        File primers = "mt_primers.fa"
        File sam = "${sample_name}.mt_primers.sam"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "40G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
    }
}
