version 1.0

workflow Bifrost_DBGWAS{
    meta{
        description: "a workflow that construct bifrost colored de bruijn graph"
    }
    input{
        String input_path
        File reference_fasta
        String outputprefix
        Int k
    }
    call extractfa{input: filepath=input_path}
    call construct{input: fas = extractfa.faFiles, ref = reference_fasta, outputpref = outputprefix, kmersize = k}
    output{
        File graph = construct.graph
        File graph_index = construct.graph_index
        File colors = construct.color_file
        File inputfasta = construct.fasta

    }
}

task extractfa{
    input{
        String filepath
    }
    command{
        gsutil ls "${filepath}" > fafiles.txt
        grep -E "\.fna$" fafiles.txt > fa_files.txt
        cat fa_files.txt
    }
    output {
        Array[String] faFiles = read_lines("fa_files.txt")
    }
    runtime {
        docker: "broadinstitute/gatk:4.4.0.0"
        disks: "local-disk 100 HDD"
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
    Bifrost build -t ~{num_threads} -k ~{kmersize} -i -d -c -s "~{sep='" -s "' fas}" -r ~{ref} -o ~{outputpref}_Bfrost_graph
    cat ~{sep=" " fas} > all.fasta
    >>>
   

    output{
        File color_file="~{outputpref}_Bfrost_graph.color.bfg"
        File graph = "~{outputpref}_Bfrost_graph.gfa.gz"
        File graph_index = "~{outputpref}_Bfrost_graph.bfi"
        File fasta = "all.fasta"
    }

    Int disk_size = 100 

    runtime {
        cpu: num_threads
        memory: "64 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/bifrost:v1"
    }
}