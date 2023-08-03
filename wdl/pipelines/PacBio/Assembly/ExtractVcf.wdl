version 1.0

workflow ExtractVcf{
    meta{
        description: "a workflow that extract vcfs in a genomic interval"
    }
    input{
        File wholegenomevcf
        File wholegenometbi
        String genomeregion
        String prefix
        Int nthreads
    }
    call extract_vcf{input: vcf_input=wholegenomevcf, vcf_index=wholegenometbi, region=genomeregion, prefix=prefix}
    output{
        File subset_vcf = extract_vcf.local_vcf
    }
}

task extract_vcf{
    input{
        File vcf_input
        File vcf_index
        String region
        String prefix
    }
    command <<<
        bcftools view -r ~{region} ~{vcf_input} -o ~{prefix}_subset.vcf
    >>>

    output{
        File local_vcf="~{prefix}_subset.vcf"
    }

    Int disk_size = 1 + ceil(2 * size(vcf_input, "GiB"))

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