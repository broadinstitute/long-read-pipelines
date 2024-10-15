version 1.0

workflow SplitbySample {
    input {
        String sample_id
        File eval_vcf
        File eval_vcf_index

    }
    call SplitVCFbySample { input:
        joint_vcf = eval_vcf,
        joint_vcf_index = eval_vcf_index,
        samplename = sample_id
    }


     output {
        File vcf = SplitVCFbySample.single_sample_vcf
        File vcf_index = SplitVCFbySample.single_sample_vcf_tbi
        
    }
}



task SplitVCFbySample {
    input{       
        File joint_vcf
        File joint_vcf_index
        String samplename
        String prefix
    }
    
    command <<<
        set -x pipefail

        time bcftools view --threads 4 -s ~{samplename} ~{joint_vcf} -o ~{prefix}.~{samplename}.subset.g.vcf.gz

        tabix -p vcf ~{prefix}.~{samplename}.subset.g.vcf.gz

    >>>
    
    output {
		File single_sample_vcf = "~{prefix}.~{samplename}.subset.g.vcf.gz"
        File single_sample_vcf_tbi = "~{prefix}.~{samplename}.subset.g.vcf.gz.tbi"
    }


    Int disk_size = 1 + ceil(2 * (size(joint_vcf, "GiB")))

    runtime {
        cpu: 4
        memory: "24 GiB"
        disks: "local-disk " + disk_size + " SSD" #"local-disk 100 HDD"
        preemptible: 1
        maxRetries: 0
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}
