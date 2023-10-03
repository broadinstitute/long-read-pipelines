version 1.0

workflow StatisticalPhasing_whole{
    meta{
        description: "a workflow that extract vcfs and bams in a genomic interval"
    }
    input{
        Array[File] wholegenomebams
        Array[File] wholegenomebai
        Array[String] sampleIds
        File joint_vcf
        File joint_vcf_tbi
        File reference
        File genetic_mapping_tsv
        String genomeregion
        String prefix
        Int num_t
    }

    Int data_length = length(wholegenomebams)
    Array[Int] indexes= range(data_length)

    scatter (idx in indexes) {
        File bam = wholegenomebams[idx]
        File bai = wholegenomebai[idx]
        String sampleid = sampleIds[idx]

        call whatshap_phasing{  input:
            inputbams = bam,
            inputbais = bai,
            ref = reference,
            joint_vcf = joint_vcf,
            joint_vcf_tbi = joint_vcf_tbi,
            region = genomeregion,
            samplename = sampleid
        }
    }

    call merge{input:
        vcf_input = whatshap_phasing.phased_vcf,
        tbi_input = whatshap_phasing.phase_vcf_tbi,
        pref = prefix

    }

    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv)

    call shapeit4 {input:
        vcf_input = merge.merged_vcf,
        vcf_index = merge.merged_tbi,
        mappingfile = genetic_mapping_dict[genomeregion],
        region = genomeregion,
        num_threads = num_t

    }

    output{
        File phased_scaffold = shapeit4.scaffold_vcf
    }
}



task whatshap_phasing {
    input{
        File inputbams
        File inputbais
        File ref
        File joint_vcf
        File joint_vcf_tbi
        String region
        String samplename
    }
    
    command <<<

        set -x pipefail

        samtools faidx ~{ref}

        bcftools view -s ~{samplename} ~{joint_vcf} -r ~{region} -o ~{samplename}.subset.g.vcf.bgz

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


task merge{
    input{
        Array[File] vcf_input
        Array[File] tbi_input
        String pref
    }
    command <<<
        bcftools merge --merge all ~{sep=" " vcf_input} -O v -o ~{pref}.AllSamples.vcf
        bgzip -c ~{pref}.AllSamples.vcf > ~{pref}.AllSamples.vcf.gz
        tabix -p vcf ~{pref}.AllSamples.vcf.gz
    >>>

    output{
        File merged_vcf = "~{pref}.AllSamples.vcf.gz"
        File merged_tbi = "~{pref}.AllSamples.vcf.gz.tbi"
    }

    Int disk_size = 10 + ceil(2 * size(vcf_input, "GiB"))

    runtime {
        cpu: 2
        memory: "8 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        # bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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