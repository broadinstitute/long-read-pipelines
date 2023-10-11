version 1.0

workflow StatisticalPhasing_phaserare{
    meta{
        description: "a workflow that using Whatshap to do read-based phasing and Shapeit for statistical phasing"
    }
    input{
        File physical_phased_vcf
        File physical_phased_tbi
        File geneticmapping
        String shapeit4genomeregion
        String shapeit5genomeregion
        Int nthreads
    }
    call shapeit4{input: vcf_input=physical_phased_vcf, vcf_index=physical_phased_tbi, mappingfile= geneticmapping,region=shapeit4genomeregion,num_threads=nthreads}
    call shapeit5_phase_rare{input: common_scaffold=shapeit4.scaffold_vcf, vcf_input=physical_phased_vcf, vcf_index=physical_phased_tbi, mappingfile= geneticmapping, genomeregion = shapeit5genomeregion, num_threads = nthreads}

    output{
        File shapeit4_scaffold = shapeit4.scaffold_vcf
        File final_scaffold = shapeit5_phase_rare.scaffold_vcf
    }
}

task shapeit4{
    input{
        File vcf_input
        File vcf_index
        File mappingfile
        String region
        Int num_threads
    }
    command <<<
    # add AN AC tag
    shapeit4 --input ~{vcf_input} --map ~{mappingfile} --region ~{region} --use-PS 0.0001 --sequencing --output scaffold.bcf --thread ~{num_threads} --log phased.log
    
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

task shapeit5_phase_rare{
    input{
        File common_scaffold
        File vcf_input
        File vcf_index
        File mappingfile
        String genomeregion
        Int num_threads
    }
    command <<<
    # add AN AC tag
    bcftools +fill-tags ~{vcf_input} -Ob -o tmp.out.bcf -- -t AN,AC
    bcftools index tmp.out.bcf
    bcftools +fill-tags ~{common_scaffold} -Ob -o tmp.scaffold.bcf -- -t AN,AC
    bcftools index tmp.scaffold.bcf
    phase_rare_static --input tmp.out.bcf --scaffold tmp.scaffold.bcf --map ~{mappingfile} --input-region ~{genomeregion} --scaffold-region ~{genomeregion} --output shapeit5_all_variants_phased.bcf --thread ~{num_threads}
    # bcftools index target.phased.bcf
    >>>

    output{
        File scaffold_vcf = "shapeit5_all_variants_phased.bcf"
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