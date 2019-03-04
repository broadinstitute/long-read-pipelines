# see https://github.com/Oshlack/STRetch/wiki/Running-STRetch

workflow StretchWorkflow {
    File reference_files_dir
    String input_bam
    String input_bam_bai
    File repeat_regions_bed
    String output_prefix  # An label that describes this run and that can be added to output filenames and/or logs.
    String stretch_image="weisburd/stretch@sha256:d9ef8dd1d13b1801483ff1204b6ee18dfff88ccc3fdd1e013bd3e83eaf6bb328"

    call Stretch {
        input:
            reference_files_dir=reference_files_dir,
            input_bam=input_bam,
            input_bam_bai=input_bam_bai,
            repeat_regions_bed=repeat_regions_bed,
            output_prefix=output_prefix,
            docker_image=stretch_image
    }
}

task Stretch {
    File reference_files_dir
    String input_bam
    String input_bam_bai
    File repeat_regions_bed
    String output_prefix
    String docker_image

    File ref_fasta = "${reference_files_dir}/full_genome_and_STRdecoys.sorted.fasta"
    File ref_fasta_fai = "${reference_files_dir}/full_genome_and_STRdecoys.sorted.fasta.fai"
    File ref_fasta_dict = "${reference_files_dir}/full_genome_and_STRdecoys.sorted.dict"
    File ref_fasta_amb = "${reference_files_dir}/full_genome_and_STRdecoys.sorted.fasta.amb"
    File ref_fasta_ann = "${reference_files_dir}/full_genome_and_STRdecoys.sorted.fasta.ann"
    File ref_fasta_bwt = "${reference_files_dir}/full_genome_and_STRdecoys.sorted.fasta.bwt"
    File ref_fasta_genome = "${reference_files_dir}/full_genome_and_STRdecoys.sorted.fasta.genome"
    File ref_fasta_pac = "${reference_files_dir}/full_genome_and_STRdecoys.sorted.fasta.pac"
    File ref_fasta_sa = "${reference_files_dir}/full_genome_and_STRdecoys.sorted.fasta.sa"
    File ref_str_decoys_bed = "${reference_files_dir}/STRdecoys.sorted.bed"

    Int cpus = 8
    Float ref_size = 15
    Int disk_size = ceil(2.1*size(input_bam, "GB") + ref_size + 10)

    command {
        ln -s ${reference_files_dir} /STRetch/reference-data
        set -xe

        free -h

        echo '
threads = 4
refdir = "${reference_files_dir}"
REF = "${reference_files_dir}/full_genome_and_STRdecoys.sorted.fasta"
DECOY_BED = "${reference_files_dir}/STRdecoys.sorted.bed"
STR_BED = "${reference_files_dir}/stretch_repeat_regions.bed"
' >> /STRetch/pipelines/pipeline_config.groovy

        python3 /compute_print_read_intervals.py --combine-intervals ${repeat_regions_bed} > print_reads_intervals.bed
        java -Xms2g -jar /gatk.jar PrintReads \
                    -R ${ref_fasta} \
                    -I ${input_bam} \
                    --interval-padding 5000 \
                    -L print_reads_intervals.bed \
                    -L UNMAPPED \
                    -O local.sharded.bam
        ln -s local.sharded.bai local.sharded.bam.bai
        bedtools sort -i ${repeat_regions_bed} > ${repeat_regions_bed}.sorted.bed
        python3 /generate_repeat_specs.py --stretch ${ref_fasta} ${repeat_regions_bed}.sorted.bed ${reference_files_dir}/stretch_repeat_regions.bed
        ln -s ${reference_files_dir}/stretch_repeat_regions.bed
        /STRetch/tools/bin/bpipe run \
            -p bwa_parallelism=2 \
            -p input_regions=${reference_files_dir}/stretch_repeat_regions.bed \
            /STRetch/pipelines/STRetch_wgs_bam_pipeline.groovy \
            local.sharded.bam
        mv STRs.tsv ${output_prefix}.STRs.tsv
    }

    output {
        File stretch_repeat_regions = "stretch_repeat_regions.bed"
        File print_reads_intervals = "print_reads_intervals.bed"
        File tmp_files = "tmp.*"
        File local_merged_bam = "local.*.merge.bam"
        File local_merged_bam_bai = "local.*.merge.*.bai"
        File locus_counts = "*.locus_counts"
        File str_counts = "*.STR_counts"
        File median_cov = "*.median_cov"
        File strs_tsv = "*.STRs.tsv"

    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "52G"
        disks: "local-disk ${disk_size} SSD"
    }
}
