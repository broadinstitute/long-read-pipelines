version 1.0

import "HybridPhaseOneChrSmallVariants.wdl" as PhaseOneChr

workflow HybridPhaseWholeGenome {
    meta{
        description: "a workflow that extract vcfs and bams in a genomic interval"
    }
    input{
        Array[File] whole_genome_bams
        Array[File] whole_genome_bais
        Array[String] sampleIds
        File joint_vcf
        File joint_vcf_tbi
        File reference
        File genetic_mapping_tsv
        String prefix
        Int num_t
    }
    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv)

    Array[String] chr_list = ["chr20", "chr21"]

    # double scatter: first by chr, then by sample
    scatter (genome_region in chr_list) {
        
        call extract_vcf { input: 
            vcf_input=joint_vcf, vcf_index=joint_vcf_tbi, region=genome_region, prefix=prefix
        }

    }

    output{
        Array[File] phased_scaffold = shapeit4.scaffold_vcf
    }
}

task extract_vcf {
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

task subset_vcf {
    input{
        File vcf_input
        File bcf_input
        String region
        String prefix
    }
    command <<<
        bcftools index ~{vcf_input}
        bcftools view -r ~{region} ~{vcf_input} -o ~{prefix}_truth_subset.vcf
        bgzip -c ~{prefix}_truth_subset.vcf > ~{prefix}_truth_subset.vcf.gz
        tabix -p vcf ~{prefix}_truth_subset.vcf.gz
        

        bcftools index ~{bcf_input}
        bcftools view -r ~{region} ~{bcf_input} -o ~{prefix}_test_subset.vcf
        bgzip -c ~{prefix}_test_subset.vcf > ~{prefix}_test_subset.vcf.gz
        tabix -p vcf ~{prefix}_test_subset.vcf.gz
    >>>

    output{
        File local_truth_vcf="~{prefix}_truth_subset.vcf.gz"
        File local_truth_tbi="~{prefix}_truth_subset.vcf.gz.tbi"
        File local_test_vcf = "~{prefix}_test_subset.vcf.gz"
        File local_test_vcf_index="~{prefix}_test_subset.vcf.gz.tbi"
    }

    Int disk_size = 100 + ceil(2 * size(vcf_input, "GiB"))

    runtime {
        cpu: 1
        memory: "10 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}
