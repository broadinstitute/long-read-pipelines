version 1.0

workflow HLAAssociation {

    input {
        File genotype_vcf
        File reference_vcf
        String prefix
        String outputdirectory
    }

    call FormatTransition { input:
        vcf = genotype_vcf,
        outputprefix = prefix
    }

    output {
        File bedfile = FormatTransition.bed
    }
}

task FormatTransition {
    input {
        File vcf
        String outputprefix
        Int memory
    }
    command <<<
        plink --vcf ~{vcf} --make-bed --out ~{outputprefix} --vcf-half-call m 
    >>>

    output {
        File bed = "~{outputprefix}.bed"
        File bim = "~{outputprefix}.bim"
        File fam = "~{outputprefix}.fam"
        File log = "~{outputprefix}.log"
        File nosex = "~{outputprefix}.nosex"
    }
    ###################
    runtime {
        cpu: 4
        memory: memory + " GiB"
        disks: "local-disk 375 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"hangsuunc/hla_tapas:v1"
    }
}

task MakeReference {
    input {
        File reference_vcf
        String outputprefix
        Int memory
    }
    command <<<
        plink --vcf ~{reference_vcf} --make-bed --out ~{outputprefix}
        plink --bfile ~{outputprefix} --hardy --allow-no-sex --out ~{outputprefix}.miss.frq.diff
        cat ~{outputprefix}.miss.frq.diff.hwe | awk '/ALL/{split($6,V,"/"); if(V[1]+V[2]+V[3]==0){Freq=0}else{Freq=(V[2]/2+V[1])/(V[1]+V[2]+V[3])}print $2, $4, $5, Freq}' > ~{outputprefix}.miss.frq.diff.hwe.freq 
        cat ~{outputprefix}.miss.frq.diff.hwe.freq | sort -k1,1 | join - <(cat ~{outputprefix}.bim | awk '{print $2, $1 "_" $4}' | sort -k1,1) | awk '{print $5, $2, $3, $4}' | sort -k1,1 > ~{outputprefix}.Ref.Frq.chr_pos_allele

    >>>

    output {
        File bed = "~{outputprefix}.bed"
        File bim = "~{outputprefix}.bim"
        File fam = "~{outputprefix}.fam"
        File log = "~{outputprefix}.log"
        File nosex = "~{outputprefix}.nosex"
        File hwe = "~{outputprefix}.miss.frq.diff.hwe"
        File freq = "~{outputprefix}.miss.frq.diff.hwe.freq"
        File pos_allele = "~{outputprefix}.Ref.Frq.chr_pos_allele"
    }
    ###################
    runtime {
        cpu: 4
        memory: memory + " GiB"
        disks: "local-disk 375 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"hangsuunc/hla_tapas:v1"
    }
}


task SNP2HLA {

    input {
        File joint_vcf
        File joint_vcf_tbi
        Int memory
    }

    command <<<
        set -euxo pipefail
        mkdir output
        bcftools +split -Oz -o output ~{joint_vcf}
        # cd output
        # for vcf in $(find . -name "*.vcf.gz"); do
        #     tabix -p vcf "$vcf"
        # done
        # cd -

    >>>

    output {
        Array[File] vcf_by_sample = glob("output/*vcf.gz")
        # Array[File] vcf_by_sample_tbi = glob("output/*vcf.gz.tbi")

    }
    ###################
    runtime {
        cpu: 4
        memory: memory + " GiB"
        disks: "local-disk 375 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"hangsuunc/hla_tapas:v1"
    }
}
