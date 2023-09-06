version 1.0

workflow MergeBam{
    meta{
        description: "a workflow that merge multiple bam in to a single one"
    }
    input{
        Array[File] input_bams
        String outputprefix
    }
    call merge{input: bam_input = input_bams, pref = outputprefix}

    output{
        File merged_bam = merge.outputbam

    }
}

task merge{
    input{
        Array[File] bam_input
        String pref
    }
    command <<<
        samtools merge -o ~{pref}.bam ~{sep=" " bam_input}
    >>>

    output{

        File outputbam="~{pref}.bam"

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
