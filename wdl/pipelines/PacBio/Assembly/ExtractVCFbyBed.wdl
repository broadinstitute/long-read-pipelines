version 1.0

workflow ExtractVcfbyBed{
    meta{
        description: "a workflow that extract vcfs in a genomic interval"
    }
    input{
        File wholegenomevcf
        File wholegenometbi
        File bedfile
        String prefix
    }
    call extract_vcf{input: vcf_input=wholegenomevcf, vcf_index=wholegenometbi, bed_input=bedfile, prefix=prefix}
    output{
        File subset_vcf = extract_vcf.local_vcf
        File subset_tbi = extract_vcf.loca_tbi
    }
}

task extract_vcf{
    input{
        File vcf_input
        File vcf_index
        File bed_input
        String prefix
    }
    command <<<
        bcftools view -R ~{bed_input} ~{vcf_input} -Oz -o ~{prefix}_subset.vcf.gz
        tabix -p vcf ~{prefix}_subset.vcf.gz
    >>>

    output{
        File local_vcf="~{prefix}_subset.vcf.gz"
        File loca_tbi="~{prefix}_subset.vcf.gz.tbi"
    }

    Int disk_size = 50 + ceil(2 * size(vcf_input, "GiB"))

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