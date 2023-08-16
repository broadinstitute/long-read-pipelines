version 1.0
import "../../../tasks/VariantCalling/CallVariantsReadBased.wdl" as VAR

workflow mapping_phasing{
    input{
        File baminput
        File reference
        String outputpref
    }
    call bcf_mpileup{input: bam_input=baminput, ref=reference, outputprefix=outputpref}
    call whatshap_phasing{input: inputbam=baminput, ref=reference, outputprefix=outputpref, vcf=bcf_mpileup.bcf_file}
    meta{
        Purpose:"Call variants via bcf mpileup, phasing based on whatshap, assemble by wtdbg"
    }
}

task bcf_mpileup {
    input{
        File bam_input
        File ref
        String outputprefix
    }
    command <<<
        bcftools mpileup -Ou -f ~{ref} ~{bam_input} | bcftools call -mv -Ob -o ~{outputprefix}.bcf
    >>>

    output {
        File bcf_file="~{outputprefix}.bcf"
    }

    Int disk_size = 1 + ceil(2 * size(bam_input, "GiB") + size(ref, "GiB"))

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

task whatshap_phasing {
    input{
        File inputbam
        File ref
        File vcf
        String outputprefix
    }
    command <<<
        samtools faidx ~{ref}
        samtools index ~{inputbam}
        bcftools index ~{vcf}
        # phasing and index vcf files
        whatshap phase -o ~{outputprefix}.phased.vcf --tag=HP --reference=~{ref} ~{vcf} ~{inputbam}
        bgzip -c ~{outputprefix}.phased.vcf > ~{outputprefix}.phased.vcf.gz
        tabix -p vcf ~{outputprefix}.phased.vcf.gz
        rm ~{outputprefix}.phased.vcf

        # tag reads
        whatshap haplotag -o ~{outputprefix}.haplotagged.bam --output-haplotag-list ~{outputprefix}.haplotaglist.tsv.gz --reference ~{ref} ~{outputprefix}.phased.vcf.gz ~{inputbam}

        #split reads
        whatshap split --output-h1 ~{outputprefix}.hap1.bam --output-h2 ~{outputprefix}.hap2.bam ~{inputbam} ~{outputprefix}.haplotaglist.tsv.gz


    >>>

    output {
        File hap1_file="~{outputprefix}.hap1.bam"
        File hap2_file="~{outputprefix}.hap2.bam"
        File taggedreads_file="~{outputprefix}.haplotagged.bam"
        File haplotaglist="~{outputprefix}.haplotaglist.tsv.gz"
    }

    Int disk_size = 1 + ceil(2 * size(inputbam, "GiB") + size(ref, "GiB"))

    runtime {
        cpu: 1
        memory: "10 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/whatshap:v1"
    }
}

task hifiasm_asm{
    input{
        File reads
        String prefix
        Int num_cpus = 32
    }

    Int disk_size = 1 + ceil(2 * size(reads, "GiB"))

    command <<<

        set -euxo pipefail

        hifiasm -o ~{prefix} -t~{num_cpus} ~{reads}
        awk '/^S/{print ">"$2; print $3}' ~{prefix}.bp.p_ctg.gfa > ~{prefix}.bp.p_ctg.fa
        awk '/^S/{print ">"$2;print $3}' ~{prefix}.bp.hap1.p_ctg.gfa > ~{prefix}.bp.hap1.p_ctg.fa
        awk '/^S/{print ">"$2;print $3}' ~{prefix}.bp.hap2.p_ctg.gfa > ~{prefix}.bp.hap2.p_ctg.fa
    >>>

    output{
        File assembly_hap1="~{prefix}.bp.hap1.p_ctg.fa"
        File assembly_hap2="~{prefix}.bp.hap2.p_ctg.fa"
    }
    runtime {
        cpu: num_cpus
        memory: "10 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/assembly:v1"
    }    
}