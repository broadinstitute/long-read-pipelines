
workflow ExpansionHunter3Workflow {
    Int num_shards
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String input_bam
    String input_bam_bai
    File repeat_regions_bed
    Int read_depth
    String output_prefix
    String expansion_hunter_image="weisburd/expansionhunter3@sha256:24825d5b302a7dd3d9994f7350ab4ecf2484d3da26a3350627734193e1716628"

    call ScatterIntervals {
        input:
            input_bed=repeat_regions_bed,
            num_shards=num_shards,
            docker_image=expansion_hunter_image
    }

    scatter (repeat_regions_bed in ScatterIntervals.output_beds) {

        call ExpansionHunter3 {
            input:
                ref_fasta=ref_fasta,
                ref_fasta_fai=ref_fasta_fai,
                ref_fasta_dict=ref_fasta_dict,
                input_bam=input_bam,
                input_bam_bai=input_bam_bai,
                repeat_regions_bed=repeat_regions_bed,
                read_depth=read_depth,
                num_shards=num_shards,
                output_prefix=output_prefix,
                docker_image=expansion_hunter_image
        }
    }

    call GatherVCFs {
        input:
            ref_fasta=ref_fasta,
            ref_fasta_fai=ref_fasta_fai,
            input_vcfs=ExpansionHunter3.output_vcf_gz,
            input_vcfs_tbi=ExpansionHunter3.output_vcf_gz_tbi,
            output_prefix=output_prefix,
            docker_image=expansion_hunter_image
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

task ExpansionHunter3 {
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String input_bam
    String input_bam_bai
    File repeat_regions_bed
    Int read_depth
    Int num_shards
    String output_prefix
    String docker_image

    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_fai, "GB")
    Int disk_size = ceil(3*size(input_bam, "GB")/num_shards + ref_size + 10)

    command {
        set -xe

        python3 /compute_print_read_intervals.py --combine-intervals ${repeat_regions_bed} > print_reads_intervals.bed
        python3 /generate_repeat_specs.py --expansion-hunter3 ${ref_fasta} ${repeat_regions_bed} variant_catalog.json

        java -Xms2g -jar /gatk.jar PrintReads \
            -R ${ref_fasta} \
            -I ${input_bam} \
            --interval-padding 5000 \
            -L print_reads_intervals.bed \
            -L UNMAPPED \
            -O local.sharded.bam

        cp local.sharded.bai local.sharded.bam.bai

        /ExpansionHunter/bin/ExpansionHunter \
            --reference ${ref_fasta} \
            --reads local.sharded.bam \
            --variant-catalog variant_catalog.json \
            --genome-coverage ${read_depth} \
            --output-prefix ${output_prefix}

        bcftools reheader --fai ${ref_fasta_fai} ${output_prefix}.vcf | bcftools sort - | bgzip -c  > ${output_prefix}.vcf.gz
        tabix ${output_prefix}.vcf.gz
    }

    output {
        File output_vcf_gz="${output_prefix}.vcf.gz"
        File output_vcf_gz_tbi="${output_prefix}.vcf.gz.tbi"
        File output_log="${output_prefix}.log"
        File output_json="${output_prefix}.json"
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

        tabix ${output_prefix}.vcf.gz

        python3 /standardize_output_vcf.py --expansion-hunter ${output_prefix}.vcf.gz ${output_prefix}.standard_format.vcf
        cat ${output_prefix}.standard_format.vcf | vt decompose -s - | vt normalize -r ${ref_fasta} - | bgzip > ${output_prefix}.normalized.vcf.gz

        python3 /update_vcf_id_field.py ${output_prefix}.normalized.vcf.gz ${output_prefix}.normalized.updated_id.vcf
        bgzip ${output_prefix}.normalized.updated_id.vcf
        tabix ${output_prefix}.normalized.updated_id.vcf.gz
    }

    output {
        File output_vcf_gz = "${output_prefix}.vcf.gz"
        File output_vcf_gz_tbi = "${output_prefix}.vcf.gz.tbi"
        File output_normalized_vcf_gz = "${output_prefix}.normalized.updated_id.vcf.gz"
        File output_normalized_vcf_gz_tbi = "${output_prefix}.normalized.updated_id.vcf.gz.tbi"
    }

    runtime {
        docker: "${docker_image}"
    }
}

