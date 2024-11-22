version 1.0

workflow Shapeit5_switch{
    meta{
        description: "a workflow that using Shapeit5 to do phasing"
    }

    input{
        File truthvcf
        File truthvcf_index
        File testvcf
        File testvcf_index
        String? extra_args # input of regions
        String output_prefix
        Int nthreads
    }
    
    call switch{input: truth_bcf=truthvcf, truth_bcf_index=truthvcf_index, test_bcf= testvcf, test_bcf_index = testvcf_index, extra_args=extra_args, outputprefix = output_prefix, num_threads=nthreads}
    
    output{
        Array[File] output_file = switch.output_files
    }
}

task switch{
    input{
        File truth_bcf
        File truth_bcf_index
        File test_bcf
        File test_bcf_index
        String outputprefix
        Int num_threads
        String? extra_args
    }
    command <<<
        switch_static \
        --validation ~{truth_bcf} \
        --estimation ~{test_bcf} \
        --output ~{outputprefix} \
        --thread ~{num_threads} \
        ~{extra_args}
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
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/shapeit5:v1"
    }
}