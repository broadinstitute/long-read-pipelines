version 1.0

workflow StatisticalPhasing_chr{
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

    call extract_vcf{input: vcf_input=joint_vcf, vcf_index=joint_vcf_tbi, region=genomeregion, prefix=prefix}

    scatter (idx in indexes) {
        File bam = wholegenomebams[idx]
        File bai = wholegenomebai[idx]
        String sampleid = sampleIds[idx]

        call whatshap_phasing{  input:
            inputbams = bam,
            inputbais = bai,
            ref = reference,
            joint_vcf = extract_vcf.local_vcf,
            joint_vcf_tbi = extract_vcf.local_tbi,
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

task extract_vcf{
    input{
        File vcf_input
        File vcf_index
        String region
        String prefix
    }
    command <<<
        bcftools view -r ~{region} ~{vcf_input} -o ~{prefix}_subset.vcf
        bgzip -c ~{prefix}_subset.vcf > ~{prefix}_subset.vcf.gz
        tabix -p vcf ~{prefix}_subset.vcf.gz
        rm ~{prefix}_subset.vcf
    >>>

    output{
        File local_vcf="~{prefix}_subset.vcf.gz"
        File local_tbi="~{prefix}_subset.vcf.gz.tbi"
    }

    Int disk_size = 1 + ceil(2 * size(vcf_input, "GiB"))

    runtime {
        cpu: 4
        memory: "16 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task extract_bam{
    input{
        File bam_input
        File bam_index
        String region
        String pref
    }
    command <<<
        samtools view --with-header ~{bam_input} -b ~{region} -o ~{pref}.bam
        samtools index ~{pref}.bam
    >>>

    output{
        File local_bam="~{pref}.bam"
        File local_bai= "~{pref}.bam.bai"
    }

    Int disk_size = 10 + ceil(2 * size(bam_input, "GiB"))

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


    Int disk_size = 100 + ceil(2 * (size(joint_vcf, "GiB")) + size(inputbams, "GiB"))

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

    Int disk_size = 100 + ceil(2 * size(vcf_input, "GiB"))

    runtime {
        cpu: 32
        memory: "128 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        # bootDiskSizeGb: 10
        preemptible: 0
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
        shapeit4 --input ~{vcf_input} --map ~{mappingfile} --region ~{region} --use-PS 0.0001 --sequencing --output ~{region}_scaffold.bcf --thread ~{num_threads} --log phased.log
    
    >>>

    output{
        File scaffold_vcf = "~{region}_scaffold.bcf"
    }

    Int disk_size = 100 + ceil(2 * size(vcf_input, "GiB"))

    runtime {
        cpu: 96
        memory: "600 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "lifebitai/shapeit4:latest"
    }
}