version 1.0

workflow StatisticalPhasing{
    meta{
        description: "a workflow that using Whatshap to do read-based phasing and Shapeit for statistical phasing"
    }
    input{
        Array[File] baminputs
        File reference
        File joint_g_vcf
        File joint_g_tbi
        String outputpref
        File geneticmapping
        String genomeregion
        Int nthreads
    }
    call whatshapphasing{input: inputbams=baminputs, ref=reference, outputprefix=outputpref, joint_vcf=joint_g_vcf, joint_vcf_tbi=joint_g_tbi}
    call shapeit5{input: vcf_input=joint_g_vcf, vcf_index=joint_g_tbi, mappingfile= geneticmapping, region=genomeregion, num_threads=nthreads}
    call shapeit4{input: vcf_input=whatshapphasing.phased_vcf, vcf_index=whatshapphasing.phase_vcf_tbi, mappingfile= geneticmapping,scaffold=shapeit5.scaffold_vcf,scaffold_index=shapeit5.scaffold_vcf_index, region=genomeregion,num_threads=nthreads}
    output{
        File scaffold = shapeit4.scaffold_vcf
    }
}

task whatshapphasing {
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
        bcftools +fill-tags ~{outputprefix}.phased.vcf -Ob -o tmp.phased.bcf -- -t AN,AC
        bcftools index tmp.phased.bcf

    >>>
    
    output {
		File phased_vcf = "tmp.phased.bcf"
        File phase_vcf_tbi = "tmp.phased.bcf.csi"
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


task shapeit5{
    input{
        File vcf_input
        File vcf_index
        File mappingfile
        String region
        Int num_threads
        Float minimal_maf = 0.001
    }
    command <<<
    # add AN AC tag
    bcftools +fill-tags ~{vcf_input} -Ob -o tmp.out.bcf -- -t AN,AC
    bcftools index tmp.out.bcf
    phase_common_static --input tmp.out.bcf --filter-maf ~{minimal_maf} --region ~{region} --map ~{mappingfile} --output scaffold.bcf --thread ~{num_threads}
    bcftools +fill-tags scaffold.bcf -Ob -o tmp.scaffold.bcf -- -t AN,AC
    bcftools index tmp.scaffold.bcf
    >>>

    output{
        File scaffold_vcf = "tmp.scaffold.bcf"
        File scaffold_vcf_index = "tmp.scaffold.bcf.csi"
    }

    Int disk_size = 100 + ceil(2 * size(vcf_input, "GiB"))

    runtime {
        cpu: 1
        memory: "50 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "lindonkambule/shapeit5_2023-05-05_d6ce1e2:v5.1.1"
    }
}

task shapeit4{
    input{
        File vcf_input
        File vcf_index
        File mappingfile
        File scaffold
        File scaffold_index
        String region
        Int num_threads
    }
    command <<<
    # add AN AC tag
    shapeit4 --input ~{vcf_input} --map ~{mappingfile} --region ~{region} --scaffold ~{scaffold} --use-PS 0.0001 --output scaffold.bcf --thread ~{num_threads} --log phased.log

    >>>

    output{
        File scaffold_vcf = "scaffold.bcf"
    }

    Int disk_size = 100 + ceil(2 * size(vcf_input, "GiB"))

    runtime {
        cpu: 1
        memory: "50 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "lifebitai/shapeit4:latest"
    }
}