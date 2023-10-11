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

task MeasurePhasingSwitchErrorRate {
    input{
        File truth_bcf
        File truth_bcf_index
        File test_bcf
        File test_bcf_index
        String region
        String outputprefix
        Int num_threads
    }
    command <<<
        switch_static --validation ~{truth_bcf} --estimation ~{test_bcf} --region ~{region} --output ~{outputprefix} --thread ~{num_threads}
    >>>

    output{
        Array[File] output_files = glob("*")
    }

    Int disk_size = 100 + ceil(2 * (size(truth_bcf, "GiB") + size(test_bcf, "GiB")))

    runtime {
        cpu: 64
        memory: "416 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "lindonkambule/shapeit5_2023-05-05_d6ce1e2:v5.1.1"
    }
}
