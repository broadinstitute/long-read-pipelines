version 1.0

workflow ConsensusAsm{
    meta{
        description: "a workflow that construct consensus assemblies from shapeit4 scaffold"
    }
    input{
        Array[String] SampleIDs
        File final_sv_scaffold
        File final_snp_scaffold
        File reference_fasta
    }

    call merge_scaffold{ input:
        input_sv_scaffold = final_sv_scaffold,
        input_snp_scaffold = final_snp_scaffold,
        prefix = "chr6"
    }

    scatter (sampleid in SampleIDs){
        call SplitVCFbySample as SP { input:
            joint_vcf = merge_scaffold.combined_scaffold,
            region = "chr6:28381449-33301940", # MHC region in T2T
            samplename = sampleid
        }
        call ConstructConsensus { input:
            single_sample_vcf = SP.single_sample_vcf,
            reference = reference_fasta,
            samplename = sampleid
        }
    }


    output{
        Array[File] Haplotype1 = ConstructConsensus.consensus_hap1
        Array[File] Haplotype2 = ConstructConsensus.consensus_hap2

    }
}

task merge_scaffold{
    input{
        File input_sv_scaffold
        File input_snp_scaffold
        File prefix
    }
    command <<<
        bcftools index ~{input_sv_scaffold}
        bcftools index ~{input_snp_scaffold}
        bcftools concat -a -D ~{input_sv_scaffold} ~{input_snp_scaffold} -O z -o ~{prefix}.combined_scaffold.bcf
        bcftools index ~{prefix}.combined_scaffold.bcf
    >>>
   

    output{
        File combined_scaffold = "~{prefix}.combined_scaffold.bcf"
        File combined_scaffold_index = "~{prefix}.combined_scaffold.bcf.csi"

    }

    Int disk_size = 100 

    runtime {
        cpu: 8
        memory: "32 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us-central1-docker.pkg.dev/broad-dsp-lrma/fusilli/fusilli:devel"
    }
}


task SplitVCFbySample {
    input{       
        File joint_vcf
        String region
        String samplename
    }
    
    command <<<
        set -x pipefail

        bcftools index ~{joint_vcf}

        bcftools view -s ~{samplename} ~{joint_vcf} -r ~{region} -o ~{samplename}.subset.g.vcf.bgz

        tabix -p vcf ~{samplename}.subset.g.vcf.bgz

    >>>
    
    output {
		File single_sample_vcf = "~{samplename}.subset.g.vcf.bgz"
        File single_sample_vcf_tbi = "~{samplename}.subset.g.vcf.bgz.tbi"
    }


    Int disk_size = 1 + ceil(2 * (size(joint_vcf, "GiB")))

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task ConstructConsensus {
    input{       
        File single_sample_vcf
        File reference
        String samplename
    }
    
    command <<<
        set -x pipefail
        bcftools index ~{single_sample_vcf}
        bcftools consensus -H 1 -f ~{reference} ~{single_sample_vcf} > ~{samplename}_MHC_hap1.fasta
        bcftools consensus -H 2 -f ~{reference} ~{single_sample_vcf} > ~{samplename}_MHC_hap2.fasta

    >>>
    
    output {
		File consensus_hap1 = "~{samplename}_MHC_hap1.fasta"
        File consensus_hap2 = "~{samplename}_MHC_hap2.fasta"
    }


    Int disk_size = 1 + ceil(2 * (size(single_sample_vcf, "GiB")))

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}