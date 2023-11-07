version 1.0
workflow DeepSomatic{
    meta{
        description: "a workflow that call somatic mutation for chrM"
    }
    input{
        File reference
        File input_baseline_bam
        File input_baseline_bam_bai
        File input_sample_bam
        File input_sample_bam_bai
        String sample_id
        String output_directory
        String chromosome
        Int nthread
    }

    call deep_somatic{ input:
        reference_fasta = reference,
        input_baseline_bam = input_baseline_bam,
        input_baseline_bam_bai = input_baseline_bam_bai,
        input_sample_bam = input_sample_bam,
        input_sample_bam_bai = input_sample_bam_bai,
        sample_id = sample_id,
        output_directory = output_directory,
        chromosome = chromosome,
        nthread = nthread
    }


    output{
        File output_vcf = deep_somatic.somatic_vcf
        File output_vcf_index = deep_somatic.somatic_vcf_index

    }
}

task deep_somatic{
    input{
        File reference_fasta
        File input_baseline_bam
        File input_baseline_bam_bai
        File input_sample_bam
        File input_sample_bam_bai
        String sample_id
        String output_directory
        String chromosome
        Int nthread
    }
    command <<<
        run_deepsomatic \
        --model_type=PACBIO \
        --ref=~{reference_fasta} \
        --reads_normal=~{input_baseline_bam} \
        --reads_tumor=~{input_sample_bam} \
        --output_vcf=~{sample_id}_deepsomatic_output.vcf.gz \
        --sample_name_tumor="baseline" \
        --sample_name_normal=~{sample_id} \
        --num_shards=~{nthread} \
        --logging_dir=~{output_directory} \
        --intermediate_results_dir=~{output_directory} \
        --regions=~{chromosome}

        bcftools index ~{sample_id}_deepsomatic_output.vcf.gz
    >>>
   

    output{
        File somatic_vcf = "~{sample_id}_deepsomatic_output.vcf.gz"
        File somatic_vcf_index = "~{sample_id}_deepsomatic_output.vcf.gz.csi"
    }

    Int disk_size = 100 

    runtime {
        cpu: 16
        memory: "64 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/deepsomatic:1.6.0"
    }
}



