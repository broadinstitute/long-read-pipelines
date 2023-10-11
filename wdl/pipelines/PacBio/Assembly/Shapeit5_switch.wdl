version 1.0

workflow Shapeit5_switch{
    meta{
        description: "a workflow that using Shapeit5 to do phasing"
    }

    input{
        File truthvcf
        File testvcf
        String genomeregion
        String output_prefix
        Int nthreads
    }
    call subset_vcf{input: vcf_input=truthvcf, bcf_input = testvcf, region = genomeregion, prefix = output_prefix}
    call switch{input: truth_bcf=subset_vcf.local_truth_vcf, truth_bcf_index=subset_vcf.local_truth_tbi, test_bcf= subset_vcf.local_test_vcf, test_bcf_index = subset_vcf.local_test_vcf_index, region=genomeregion, outputprefix = output_prefix, num_threads=nthreads}
    
    output{
        Array[File] output_file = switch.output_files
    }
}

task subset_vcf{
    input{
        File vcf_input
        File bcf_input
        String region
        String prefix
    }
    command <<<
        bcftools index ~{vcf_input}
        bcftools view -r ~{region} ~{vcf_input} -o ~{prefix}_truth_subset.vcf
        bgzip -c ~{prefix}_truth_subset.vcf > ~{prefix}_truth_subset.vcf.gz
        tabix -p vcf ~{prefix}_truth_subset.vcf.gz
        

        bcftools index ~{bcf_input}
        bcftools view -r ~{region} ~{bcf_input} -o ~{prefix}_test_subset.vcf
        bgzip -c ~{prefix}_test_subset.vcf > ~{prefix}_test_subset.vcf.gz
        tabix -p vcf ~{prefix}_test_subset.vcf.gz
    >>>

    output{
        File local_truth_vcf="~{prefix}_truth_subset.vcf.gz"
        File local_truth_tbi="~{prefix}_truth_subset.vcf.gz.tbi"
        File local_test_vcf = "~{prefix}_test_subset.vcf.gz"
        File local_test_vcf_index="~{prefix}_test_subset.vcf.gz.tbi"
    }

    Int disk_size = 100 + ceil(2 * size(vcf_input, "GiB"))

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


task switch{
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
        preemptible: 2
        maxRetries: 1
        docker: "lindonkambule/shapeit5_2023-05-05_d6ce1e2:v5.1.1"
    }
}