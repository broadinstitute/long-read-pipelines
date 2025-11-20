version 1.0

workflow IGVScreenshotWorkflow {
    
    input {
        File aligned_bam_hap1
        File aligned_bam_hap1_bai
        File aligned_bam_hap2
        File aligned_bam_hap2_bai
        File alignments
        File alignments_bai
        File bed_file
        File fasta_file
        File fasta_file_fai
        String sample_name
        Int image_height = 500
        Int memory_mb = 4000
        Int disk_gb = 100           # Disk size in GB, default to 100 GB
        String docker_image = "us.gcr.io/broad-dsp-lrma/igv_screenshot_docker:v982024"  # The Docker image to use
    }

    call RunIGVScreenshot {
        input:
            aligned_bam_hap1 = aligned_bam_hap1,
            aligned_bam_hap1_bai = aligned_bam_hap1_bai,
            aligned_bam_hap2 = aligned_bam_hap2,
            aligned_bam_hap2_bai = aligned_bam_hap2_bai,
            alignments = alignments,
            alignments_bai = alignments_bai,
            bed_file = bed_file,
            fasta_file = fasta_file,
            fasta_file_fai = fasta_file_fai,
            sample_name = sample_name,
            image_height = image_height,
            memory_mb = memory_mb,
            disk_gb = disk_gb,
            docker_image = docker_image
    }

    output {
        Array[File] snapshots = RunIGVScreenshot.snapshots
    }
}

task RunIGVScreenshot {
    
    input {
        File aligned_bam_hap1
        File aligned_bam_hap1_bai
        File aligned_bam_hap2
        File aligned_bam_hap2_bai
        File alignments
        File alignments_bai
        File bed_file
        File fasta_file
        File fasta_file_fai
        String sample_name
        Int image_height
        Int memory_mb
        Int disk_gb
        String docker_image
    }

    command <<<
        set -euo pipefail

        # Ensure the snapshots directory exists
        mkdir -p 'output/IGV_Snapshots'

        # Start a virtual frame buffer to allow IGV to render
        Xvfb :1 -screen 0 1024x768x16 &> xvfb.log &
        export DISPLAY=:1

        # Run the IGV screenshot script with the provided inputs
        python3 /opt/IGV_Linux_2.18.2/make_igv_screenshot.py \
          ~{aligned_bam_hap1} ~{aligned_bam_hap2} ~{alignments} \
          -r ~{bed_file} \
          -ht ~{image_height} \
          -bin /opt/IGV_Linux_2.18.2/igv.sh \
          -mem ~{memory_mb} \
          --fasta_file ~{fasta_file} \
          --sample_name ~{sample_name}

        # Move the screenshots to the IGV_Snapshots directory
        #mv -- *.png 'output/IGV_Snapshots/'
    >>>

    runtime {
        docker: docker_image
        memory: "~{memory_mb} MB"
        cpu: 2
        disks: "local-disk ~{disk_gb} SSD"
    }

    output {
        Array[File] snapshots = glob("output/IGV_Snapshots/*.png")
    }
}
