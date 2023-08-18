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
        File readphased_vcf_tbi = phasing.phase_vcf_tbi
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
        
        whatshap phase -o ~{outputprefix}.phased.vcf --tag=PS --reference=~{ref} ~{joint_vcf} ~{sep=" " inputbams}

        bgzip -c ~{outputprefix}.phased.vcf > ~{outputprefix}.phased.vcf.gz

        tabix -p vcf ~{outputprefix}.phased.vcf.gz

    >>>
    
    output {
		File phased_vcf = "~{outputprefix}.phased.vcf.gz"
        File phase_vcf_tbi = "~{outputprefix}.phased.vcf.gz.tbi"
    }


    Int disk_size = 100

    runtime {
        cpu: 16
        memory: "64 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/whatshap:v1"
    }
}
