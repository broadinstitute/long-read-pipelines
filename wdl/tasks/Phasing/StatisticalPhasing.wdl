version 1.0

task Shapeit4 {
    input{
        File vcf_input
        File vcf_index
        File mappingfile
        String region
        Int num_threads
    }
    command <<<
        # add AN AC tag
        shapeit4 --input ~{vcf_input} --map ~{mappingfile} --region ~{region} --use-PS 0.0001 --sequencing --output ~{region}_scaffold.bcf --thread ~{num_threads} --log phased.log
    
    >>>

    output{
        File scaffold_vcf = "~{region}_scaffold.bcf"
    }

    Int disk_size = 100 + ceil(2 * size(vcf_input, "GiB"))

    runtime {
        cpu: 96
        memory: "600 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "lifebitai/shapeit4:latest"
    }
}
