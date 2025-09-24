version 1.0

workflow IGVScreenshotWorkflow {

    input {
        File aligned_bam1
        File aligned_bam1_bai
        File? aligned_bam2
        File? aligned_bam2_bai
        File regions_bed
        File reference_fasta
        File reference_fasta_fai
        File? truth_haplotype_1
        File? truth_haplotype_1_bai
        File? truth_haplotype_2
        File? truth_haplotype_2_bai
        File? haplotype_8x_hap1
        File? haplotype_8x_hap1_bai
        File? haplotype_8x_hap2
        File? haplotype_8x_hap2_bai
        File? TRGT_VCF
        File? TRGT_VCF_tbi
        String sample_name
        Int image_height = 1000
        Int memory_mb = 4000
        Int disk_gb = 100
        String docker_image = "us.gcr.io/broad-dsp-lrma/igv_screenshot_docker:v9172024"
    }

    call RunIGVScreenshot {
        input:
            aligned_bam1 = aligned_bam1,
            aligned_bam1_bai = aligned_bam1_bai,
            aligned_bam2 = aligned_bam2,
            aligned_bam2_bai = aligned_bam2_bai,
            regions_bed = regions_bed,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            truth_haplotype_1 = truth_haplotype_1,
            truth_haplotype_1_bai = truth_haplotype_1_bai,
            truth_haplotype_2 = truth_haplotype_2,
            truth_haplotype_2_bai = truth_haplotype_2_bai,
            haplotype_8x_hap1 = haplotype_8x_hap1,
            haplotype_8x_hap1_bai = haplotype_8x_hap1_bai,
            haplotype_8x_hap2 = haplotype_8x_hap2,
            haplotype_8x_hap2_bai = haplotype_8x_hap2_bai,
            TRGT_VCF = TRGT_VCF,
            TRGT_VCF_tbi = TRGT_VCF_tbi,
            sample_name = sample_name,
            image_height = image_height,
            memory_mb = memory_mb,
            disk_gb = disk_gb,
            docker_image = docker_image
    }

    output {
        Array[File] screenshots = RunIGVScreenshot.screenshots
    }
}

task RunIGVScreenshot {

    input {
        File aligned_bam1
        File aligned_bam1_bai
        File? aligned_bam2
        File? aligned_bam2_bai
        File regions_bed
        File reference_fasta
        File reference_fasta_fai
        File? truth_haplotype_1
        File? truth_haplotype_1_bai
        File? truth_haplotype_2
        File? truth_haplotype_2_bai
        File? haplotype_8x_hap1
        File? haplotype_8x_hap1_bai
        File? haplotype_8x_hap2
        File? haplotype_8x_hap2_bai
        File? TRGT_VCF
        File? TRGT_VCF_tbi
        String sample_name
        Int image_height
        Int memory_mb
        Int disk_gb
        String docker_image
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
            ${aligned_bam1} \
            ~{if defined(aligned_bam2) then "--second_alignment_reads " + aligned_bam2 else ""} \
            ~{if defined(truth_haplotype_1) then "--truth_haplotype_1 " + truth_haplotype_1 else ""} \
            ~{if defined(truth_haplotype_2) then "--truth_haplotype_2 " + truth_haplotype_2 else ""} \
            ~{if defined(TRGT_VCF) then "--targeted_vcf " + TRGT_VCF else ""} \
            ~{if defined(haplotype_8x_hap1) then "--second_alignment_reads " + haplotype_8x_hap1 else ""} \
            ~{if defined(haplotype_8x_hap2) then "--second_alignment_reads " + haplotype_8x_hap2 else ""} \
            -r ${regions_bed} \
            -f ${reference_fasta} \
            --sample_name ${sample_name} \
            --snapshot_format png \
            --output_dir igv_output \
            -ht ${image_height}
    >>>

    runtime {
        docker: docker_image
        memory: "${memory_mb} MB"
        cpu: 2
        disks: "local-disk ${disk_gb} SSD"
    }

    output {
        Array[File] screenshots = glob("igv_output/*.png")
    }
}
