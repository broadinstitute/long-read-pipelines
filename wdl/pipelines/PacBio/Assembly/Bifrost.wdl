version 1.0

workflow Bifrost{
    meta{
        description: "a workflow that construct bifrost colored de bruijn graph"
    }
    input{
        Array[File] input_fastas
        File reference_fasta
        String outputprefix
        Int k
    }
    call construct{input: fas = input_fastas, ref = reference_fasta, outputpref = outputprefix, kmersize = k}
    output{
        File graph = construct.graph
        File colors = construct.color_file
    }
}

task construct{
    input{
        Array[File] fas
        File ref
        String outputpref
        Int num_threads = 16
        Int kmersize
    }
    command <<<
    set -x pipefail
    Bifrost build -t ~{num_threads} -k ~{kmersize} -i -d -c -s ~{sep=" " fas} -r ~{ref} -o ~{outputpref}_Bfrost_graph
    >>>

    output{
        File color_file="~{outputpref}_Bfrost_graph.bfg_colors"
        File graph = "~{outputpref}_Bfrost_graph.gfa"
    }

    Int disk_size = 100 

    runtime {
        cpu: num_threads
        memory: "10 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/bifrost:v1"
    }
}
