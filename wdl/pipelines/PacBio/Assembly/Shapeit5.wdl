version 1.0

workflow Shapeit5{
    meta{
        description: "a workflow that using Shapeit5 to do phasing"
    }
    input{
        File wholegenomevcf
        File wholegenometbi
        File geneticmapping
        String genomeregion
        Int nthreads
    }
    call phasing{input: vcf_input=wholegenomevcf, vcf_index=wholegenometbi, mappingfile= geneticmapping, region=genomeregion, num_threads=nthreads}
    output{
        File scaffold = phasing.scaffold_vcf
    }
}

task phasing{
    input{
        File vcf_input
        File vcf_index
        File mappingfile
        String region
        Int num_threads
        Float minimal_maf = 0.001
    }
    command <<<
    #bcftools +fill-tags bcftools view -Oz -o myVCF.filtag.vcf.gz
    # add AN AC tag
    bcftools +fill-tags ~{vcf_input} -Ob -o tmp.out.bcf -- -t AN,AC
    bcftools index tmp.out.bcf
    phase_common_static --input tmp.out.bcf --filter-maf ~{minimal_maf} --region ~{region} --map ~{mappingfile} --output scaffold.bcf --thread ~{num_threads}
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
        docker: "lindonkambule/shapeit5_2023-05-05_d6ce1e2:v5.1.1"
    }
}