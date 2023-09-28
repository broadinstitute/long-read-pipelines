version 1.0

workflow MergeVCFs{
    meta{
        description: "a workflow to merge vcfs"
    }
    input{
        Array[File] vcfs
        Array[File] tbis
        String Pref
        
    }
    call merge{input: vcf_input=vcfs, tbi_input=tbis, pref=Pref}
    output{
        File vcf = merge.merged_vcf
        File tbi = merge.merged_tbi
    }
}

task merge{
    input{
        Array[File] vcf_input
        Array[File] tbi_input
        String pref
    }
    command <<<
        bcftools merge --merge all ~{sep=" " vcf_input} -O v -o ~{pref}.AllSamples.vcf
        bgzip -c ~{pref}.AllSamples.vcf > ~{pref}.AllSamples.vcf.gz
        tabix -p vcf ~{pref}.AllSamples.vcf.gz
    >>>

    output{
        File merged_vcf = "~{pref}.AllSamples.vcf.gz"
        File merged_tbi = "~{pref}.AllSamples.vcf.gz.tbi"
    }

    Int disk_size = 10 + ceil(2 * size(vcf_input, "GiB"))

    runtime {
        cpu: 2
        memory: "8 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        # bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}


