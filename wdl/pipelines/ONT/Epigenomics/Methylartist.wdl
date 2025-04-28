version 1.0

workflow methylartist{
    meta{
        description: "a workflow that using methylartist to plot methylated signal from MM/ML tagged bams"
    }

    input{
        Array[File] bams
        Array[File] bais
        File reference
        String region
        String motif
    }
    call Methylartist{input: bams=bams, bais=bais, reference = reference, region = region, motif = motif}
    
    output{
        Array[File] output_file = Methylartist.output_files
    }
}

task Methylartist{
    input{
        Array[File] bams
        Array[File] bais
        File reference
        String region
        String motif = "CG"
    }

    command <<<
        methylartist locus -b ~{sep="," bams} -i ~{region} --ref ~{reference} --motif ~{motif}
    >>>

    output{
        Array[File] output_files = glob("*.png")
    }

    Int disk_size = 100 + ceil(2 * (size(bams, "GiB")))

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/methylartist:v1"
    }
}