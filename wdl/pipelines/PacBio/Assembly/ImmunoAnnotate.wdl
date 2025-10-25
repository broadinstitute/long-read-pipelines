version 1.0

workflow ImmunoAnnotate{
    meta{
        description: "a workflow that extract vcfs in a genomic interval"
    }
    input{
        File target_assembly_Hap1
        File target_assembly_Hap2
        String prefix
    }
    call immunoAnnotate as hap1_annotation {input: target_asm=target_assembly_Hap1, prefix = prefix + "hap1"}
    call immunoAnnotate as hap2_annotation {input: target_asm=target_assembly_Hap2, prefix = prefix + "hap2"}
    
    output {
        File annotation1 = hap1_annotation.gtf_file
        File annotation2 = hap2_annotation.gtf_file
    }
}

task immunoAnnotate{
    input{
        File target_asm
        String prefix
        Int threads = 2
    }
    command <<<
        bash /Immuannot/scripts.pub.v3/immuannot.sh -c ~{target_asm} -r /Immuannot/Data-2024Feb02 -o ~{prefix} -t ~{threads}
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