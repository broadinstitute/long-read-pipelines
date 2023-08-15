version 1.0

workflow WhatshapPhasing{
    input{
        Array[File] baminputs
        File reference
        File joint_g_vcf
        File joint_g_tbi
        String outputpref
    }
    call phasing{input: inputbams=baminputs, ref=reference, outputprefix=outputpref, joint_vcf=joint_g_vcf, joint_vcf_tbi=joint_g_tbi}
    meta{
        Purpose:" read-based phasing by whatshap"
    }

    output{
        File readphased_vcf = phasing.phased_vcf
    }
}


task phasing {
    input{
        Array[File] inputbams
        File ref
        File joint_vcf
        File joint_vcf_tbi
        String outputprefix
    }
    command <<<

        set -x pipefail
        samtools faidx ~{ref}
        samtools index -M ~{sep=' ' inputbams}        
        # phasing and index vcf files
        #echo "~{sep="\n" inputbams}" > list.txt
        # whatshap phase -o ~{outputprefix}.phased.vcf --tag=HP --reference=~{ref} ~{joint_vcf} ~{sep=" " inputbams}
        # bgzip -c ~{outputprefix}.phased.vcf > ~{outputprefix}.phased.vcf.gz
        # tabix -p vcf ~{outputprefix}.phased.vcf.gz
        # rm ~{outputprefix}.phased.vcf

        # tag reads
        #whatshap haplotag -o ~{outputprefix}.haplotagged.bam --output-haplotag-list ~{outputprefix}.haplotaglist.tsv.gz --reference ~{ref} ~{outputprefix}.phased.vcf.gz ~{inputbams}

        #split reads
        #whatshap split --output-h1 ~{outputprefix}.hap1.bam --output-h2 ~{outputprefix}.hap2.bam ~{inputbams} ~{outputprefix}.haplotaglist.tsv.gz


    >>>

    output {
        File phased_vcf = "~{outputprefix}.phased.vcf.gz"
        #File hap1_file="~{outputprefix}.hap1.bam"
        #File hap2_file="~{outputprefix}.hap2.bam"
        #File taggedreads_file="~{outputprefix}.haplotagged.bam"
        #File haplotaglist="~{outputprefix}.haplotaglist.tsv.gz"
    }

    Int disk_size = 100 + ceil(5 * size(inputbams, "GiB") + size(ref, "GiB"))

    runtime {
        cpu: 16
        memory: "100 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/whatshap:v1"
    }
}

