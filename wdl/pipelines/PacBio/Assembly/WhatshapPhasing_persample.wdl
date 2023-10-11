version 1.0

workflow WhatshapPhasing{
    input{
        File baminputs
        File reference
        File joint_g_vcf
        File joint_g_tbi
        String outputpref
    }
    
    
    call phasing{input: inputbams=baminputs, ref=reference, joint_vcf=joint_g_vcf, joint_vcf_tbi=joint_g_tbi, samplename=outputpref}
    
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
        File inputbams
        File ref
        File joint_vcf
        File joint_vcf_tbi
        String samplename
    }
    
    command <<<

        set -x pipefail
        samtools faidx ~{ref}
        
        samtools index -M ~{inputbams}   

        bcftools view -s ~{samplename} ~{joint_vcf} -o ~{samplename}.subset.g.vcf.bgz

        tabix -p vcf ~{samplename}.subset.g.vcf.bgz
        
        whatshap phase -o ~{samplename}.phased.vcf --tag=PS --reference=~{ref} ~{samplename}.subset.g.vcf.bgz ~{inputbams}

        bgzip -c ~{samplename}.phased.vcf > ~{samplename}.phased.vcf.gz

        tabix -p vcf ~{samplename}.phased.vcf.gz

    >>>
    
    output {
		File phased_vcf = "~{samplename}.phased.vcf.gz"
        File phase_vcf_tbi = "~{samplename}.phased.vcf.gz.tbi"
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
