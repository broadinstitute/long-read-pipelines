version 1.0

workflow IGVScreenshotWorkflow {

    input {
        File bam_file
        File bam_file_bai
        File regions_bed
        File reference_fasta
        File reference_fasta_fai
        String sample_name
        Int image_height = 1000
        Int memory_mb = 4000
        Int disk_gb = 100
        String docker_image = "us.gcr.io/broad-dsp-lrma/igv_screenshot_docker:v9172024"
        File? truth_haplotype_1
        File? truth_haplotype_1_bai
        File? truth_haplotype_2
        File? truth_haplotype_2_bai
        File? targeted_vcf
        File? targeted_vcf_tbi
        File? second_alignment_reads
        File? second_alignment_reads_bai
    }

    call RunIGVScreenshot {
        input:
            bam_file=bam_file,
            bam_file_bai=bam_file_bai,
            regions_bed=regions_bed,
            reference_fasta=reference_fasta,
            reference_fasta_fai=reference_fasta_fai,
            sample_name=sample_name,
            image_height=image_height,
            memory_mb=memory_mb,
            disk_gb=disk_gb,
            docker_image=docker_image,
            truth_haplotype_1=truth_haplotype_1,
            truth_haplotype_1_bai=truth_haplotype_1_bai,
            truth_haplotype_2=truth_haplotype_2,
            truth_haplotype_2_bai=truth_haplotype_2_bai,
            targeted_vcf=targeted_vcf,
            targeted_vcf_tbi=targeted_vcf_tbi,
            second_alignment_reads=second_alignment_reads,
            second_alignment_reads_bai=second_alignment_reads_bai
    }

    output {
        File igv_output_zip = RunIGVScreenshot.igv_output_zip
    }
}

task RunIGVScreenshot {
    
    input {
        File bam_file
        File bam_file_bai
        File regions_bed
        File reference_fasta
        File reference_fasta_fai
        String sample_name
        Int image_height
        Int memory_mb
        Int disk_gb
        String docker_image
        File? truth_haplotype_1
        File? truth_haplotype_1_bai
        File? truth_haplotype_2
        File? truth_haplotype_2_bai
        File? targeted_vcf
        File? targeted_vcf_tbi
        File? second_alignment_reads
        File? second_alignment_reads_bai
    }

    command <<<
        set -euo pipefail

        # Ensure the output directory exists
        mkdir -p igv_output

        # Start a virtual frame buffer to allow IGV to render
        Xvfb :1 -screen 0 1024x768x16 &> xvfb.log &
        export DISPLAY=:1

        # Run the IGV screenshot script with the provided inputs
        python3 /opt/IGV_Linux_2.18.2/make_igv_screenshot.py \
            ${bam_file} \
            -r ${regions_bed} \
            -f ${reference_fasta} \
            --sample_name ${sample_name} \
            --snapshot_format png \
            --output_dir igv_output \
            -ht ${image_height} \
            ~{if defined(truth_haplotype_1) then "--truth_haplotype_1 " + truth_haplotype_1 else ""} \
            ~{if defined(truth_haplotype_2) then "--truth_haplotype_2 " + truth_haplotype_2 else ""} \
            ~{if defined(targeted_vcf) then "--targeted_vcf " + targeted_vcf else ""} \
            ~{if defined(second_alignment_reads) then "--second_alignment_reads " + second_alignment_reads else ""}
        
        # Zip the output directory
        zip -r igv_output.zip igv_output/
    >>>

    runtime {
        docker: docker_image
        memory: "${memory_mb} MB"
        cpu: 2
        disks: "local-disk ${disk_gb} HDD"
    }

    output {
        File igv_output_zip = "igv_output.zip"
    }
}
