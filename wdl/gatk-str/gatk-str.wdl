
workflow GatkStrWorkflow {
    Int num_shards
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String input_bam
    String input_bam_bai
    File repeat_regions_bed
    Boolean use_str_model=true
    String output_prefix
    String gatk_str_image="weisburd/gatk-str@sha256:0708f8e42781e178503faab5f6ab7c134745f6e5bf5775d64455c53c53feabbd"

    call ScatterIntervals {
        input:
            input_bed=repeat_regions_bed,
            num_shards=num_shards,
            docker_image=gatk_str_image
    }

    scatter (repeat_regions_bed in ScatterIntervals.output_beds) {

        call GatkStr {
            input:
                ref_fasta=ref_fasta,
                ref_fasta_fai=ref_fasta_fai,
                ref_fasta_dict=ref_fasta_dict,
                input_bam=input_bam,
                input_bam_bai=input_bam_bai,
                repeat_regions_bed=repeat_regions_bed,
                num_shards=num_shards,
                use_str_model=use_str_model,
                output_prefix=output_prefix,
                docker_image=gatk_str_image
        }
    }

    call GatherVCFs {
        input:
            ref_fasta=ref_fasta,
            ref_fasta_fai=ref_fasta_fai,
            input_vcfs=GatkStr.output_normalized_vcf_gz,
            input_vcfs_tbi=GatkStr.output_normalized_vcf_gz_tbi,
            output_prefix=output_prefix,
            docker_image=gatk_str_image
    }
}

task ScatterIntervals {
    File input_bed
    Int num_shards
    String docker_image

    command {
        python3 /scatter_intervals_by_covered_region_size.py --num-shards ${num_shards} --output-dir-suffix "_sharded" ${input_bed}
    }

    output {
        Array[File] output_beds = glob("*_sharded/*.bed")
    }

    runtime {
        docker: "${docker_image}"
    }
}

task GatkStr {
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String input_bam
    String input_bam_bai
    File repeat_regions_bed
    Boolean use_str_model
    Int num_shards
    String output_prefix
    String docker_image

    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_fai, "GB")
    Int disk_size = ceil(3*size(input_bam, "GB")/num_shards + ref_size + 10)

    command {
        set -xe

        python3 /compute_print_read_intervals.py --combine-intervals ${repeat_regions_bed} > print_reads_intervals.bed
        java -Xms2g -jar /gatk.jar PrintReads \
            -R ${ref_fasta} \
            -I ${input_bam} \
            --interval-padding 5000 \
            -L print_reads_intervals.bed \
            -L UNMAPPED \
            -O local.sharded.bam
        ln -s local.sharded.bai local.sharded.bam.bai
        java -jar ${true="/GATK_withSTRmodel.jar -T HaplotypeCaller --strModelFile /wgs-all.tab" false="/gatk.jar HaplotypeCaller" use_str_model} \
            -R ${ref_fasta} \
            -I local.sharded.bam \
            -L ${repeat_regions_bed} \
            ${true="--interval_padding" false="--interval-padding" use_str_model} 300 \
            ${true="-o" false="--output" use_str_model} ${output_prefix}.vcf.gz
        tabix -f ${output_prefix}.vcf.gz
        zcat ${output_prefix}.vcf.gz | vt decompose -s - | vt normalize -r ${ref_fasta} - | bgzip > ${output_prefix}.normalized.vcf.gz
        tabix -f ${output_prefix}.normalized.vcf.gz
    }

    # && python3 /update_vcf_id_field.py ${output_prefix}.normalized.vcf.gz ${output_prefix}.normalized.updated_id.vcf

    output {
        File print_reads_intervals = "print_reads_intervals.bed"
        File output_vcf_gz="${output_prefix}.vcf.gz"
        File output_vcf_gz_tbi="${output_prefix}.vcf.gz.tbi"
        File output_normalized_vcf_gz = "${output_prefix}.normalized.vcf.gz"
        File output_normalized_vcf_gz_tbi = "${output_prefix}.normalized.vcf.gz.tbi"
    }

     runtime {
         docker: "${docker_image}"
         memory: "3 GB"
         disks: "local-disk ${disk_size} SSD"
     }
}

task GatherVCFs {
    File ref_fasta
    File ref_fasta_fai
    Array[File] input_vcfs
    Array[File] input_vcfs_tbi
    String output_prefix
    String docker_image

    command {
        bcftools concat --allow-overlaps --remove-duplicates ${sep=' ' input_vcfs} | bgzip -c > ${output_prefix}.vcf.gz
        tabix -f ${output_prefix}.vcf.gz
        zcat ${output_prefix}.vcf.gz | vt decompose -s - | vt normalize -r ${ref_fasta} - | bgzip > ${output_prefix}.normalized.vcf.gz
        tabix -f ${output_prefix}.normalized.vcf.gz

    }

    output {
        File output_vcf_gz = "${output_prefix}.vcf.gz"
        File output_vcf_gz_tbi = "${output_prefix}.vcf.gz.tbi"
        File output_normalized_vcf_gz="${output_prefix}.normalized.vcf.gz"
        File output_normalized_vcf_gz_tbi="${output_prefix}.normalized.vcf.gz.tbi"
    }

    runtime {
        docker: "${docker_image}"
    }
}
