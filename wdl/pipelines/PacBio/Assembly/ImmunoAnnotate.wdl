version 1.0

workflow ImmunoAnnotate{
    meta{
        description: "a workflow that extract vcfs in a genomic interval"
    }
    input{
        File target_asm
        String prefix
    }
    call immunoAnnotate as annotation {input: target_asm=target_asm, prefix = prefix}
    
    output {
        File output_file = annotation.gtf_file
    }
}

task immunoAnnotate{
    input{
        File target_asm
        String prefix
        Int threads = 2
    }
    command <<<
        gunzip -c ~{target_asm} > target_asm.fasta
        bash /Immuannot/scripts.pub.v3/immuannot.sh -c target_asm.fasta -r /Immuannot/Data-2024Feb02 -o ~{prefix} -t ~{threads}
    >>>
    
    output{
        File gtf_file="~{prefix}.gtf.gz"
        
    }

    Int disk_size = 1 + ceil(2 * size(target_asm, "GiB"))

    runtime {
        cpu: 2
        memory: "8 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        # bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/immuannot:v1"
    }
}