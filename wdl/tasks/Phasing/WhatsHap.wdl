version 1.0

task Phase {
    input{
        File bam
        File bai

        File ref
        
        File joint_vcf
        File joint_vcf_tbi
        
        String region
        String samplename
    }
    
    command <<<
        set -x pipefail

        samtools faidx ~{ref}

        time \
        bcftools view -s ~{samplename} ~{joint_vcf} -r ~{region} -o ~{samplename}.subset.g.vcf.bgz

        tabix -p vcf ~{samplename}.subset.g.vcf.bgz
        
        whatshap phase -o ~{samplename}.phased.vcf --tag=PS --reference=~{ref} ~{samplename}.subset.g.vcf.bgz ~{bam}

        bgzip -c ~{samplename}.phased.vcf > ~{samplename}.phased.vcf.gz

        tabix -p vcf ~{samplename}.phased.vcf.gz

    >>>
    
    output {
		File phased_vcf = "~{samplename}.phased.vcf.gz"
        File phase_vcf_tbi = "~{samplename}.phased.vcf.gz.tbi"
    }


    Int disk_size = 100 + ceil(2 * (size(joint_vcf, "GiB")) + size(bam, "GiB"))

    runtime {
        cpu: 16
        memory: "64 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "hangsuunc/whatshap:v1"
    }
}
