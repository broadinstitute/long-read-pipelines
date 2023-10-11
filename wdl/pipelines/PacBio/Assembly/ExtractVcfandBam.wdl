version 1.0

workflow ExtractVcfandBam{
    meta{
        description: "a workflow that extract vcfs and bams in a genomic interval"
    }
    input{
        File wholegenomebam
        File wholegenomebai
        File wholegenomevcf
        File wholegenometbi
        String genomeregion
        String prefix
    }
    call extract_vcf{input: vcf_input=wholegenomevcf, vcf_index=wholegenometbi, region=genomeregion, prefix=prefix}
    call extract_bam{input: bam_input=wholegenomebam, bam_index=wholegenomebai, region= genomeregion, pref=prefix}
    output{
        File subset_vcf = extract_vcf.local_vcf
        File subset_tbi = extract_vcf.local_tbi
        File subset_bam = extract_bam.local_bam
        String read_number = read_string(extract_bam.read_number) 
        #Map[String, Pair[File, File]] data = {prefix:(subset_vcf, subset_bam)}
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
        bgzip -c ~{prefix}_subset.vcf > ~{prefix}_subset.vcf.gz
        tabix -p vcf ~{prefix}_subset.vcf.gz
        rm ~{prefix}_subset.vcf
    >>>

    output{
        File local_vcf="~{prefix}_subset.vcf.gz"
        File local_tbi="~{prefix}_subset.vcf.gz.tbi"
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

task extract_bam{
    input{
        File bam_input
        File bam_index
        String region
        String pref
    }
    command <<<
        samtools view --with-header ~{bam_input} -b ~{region} -o ~{pref}.bam
        samtools view -c ~{pref}.bam > read_num.txt
        #samtools index ~{pref}.~{region}.bam
    >>>

    output{
        File local_bam="~{pref}.bam"
        File read_number = "read_num.txt"
        #File local_bai="~{pref}.~{region}.bai"
    }

    Int disk_size = 10 + ceil(2 * size(bam_input, "GiB"))

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
