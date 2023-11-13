version 1.0

workflow BandagePlotGraph{
    meta{
        description: "a workflow that using Bandage to plot gfa graph"
    }
    input{
        File graph_gfa
        String prefix
    }
    call bandage{input: gfa=graph_gfa, outputprefix=prefix}
    output{
        File figure = bandage.image
    }
}

task bandage{
    input{
        File gfa
        String outputprefix
    }
    command <<<
        Bandage image ~{gfa} ~{outputprefix}.png
    >>>

    output{
        File image = "~{outputprefix}.png"
    }

    Int disk_size = 10 + ceil(size(gfa, "GiB"))

    runtime {
        cpu: 4
        memory: "16 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "hangsuunc/bandage:v1"
    }
}