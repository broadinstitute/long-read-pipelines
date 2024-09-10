version 1.0

workflow IGVScreenshotWorkflow {

    input {
        File aligned_bam_hap1       # BAM file for haplotype 1
        File aligned_bam_hap2       # BAM file for haplotype 2
        File alignments             # BAM file for total alignments
        File bed_file               # BED file with regions
        File fasta_file             # Reference FASTA file
        File fasta_fai              # FASTA index (.fai) file
        String sample_name          # Sample name to use in filenames
        Int image_height = 500
        Int memory_mb = 4000
        Int disk_gb = 100           # Disk size in GB, default to 100 GB
        String docker_image = "us.gcr.io/broad-dsp-lrma/igv_screenshot_docker:v982024"  # The Docker image to use
    }

    call RunIGVScreenshot {
        input:
            aligned_bam_hap1 = aligned_bam_hap1,
            aligned_bam_hap1_bai = aligned_bam_hap1 + ".bai",  # Automatically infer BAI location
            aligned_bam_hap2 = aligned_bam_hap2,
            aligned_bam_hap2_bai = aligned_bam_hap2 + ".bai",  # Automatically infer BAI location
            alignments = alignments,
            alignments_bai = alignments + ".bai",  # Automatically infer BAI location
            bed_file = bed_file,
            fasta_file = fasta_file,
            fasta_fai = fasta_fai,
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
        File fasta_fai
        String sample_name
        Int image_height
        Int memory_mb
        Int disk_gb
        String docker_image
    }

    command {
        mkdir -p IGV_Snapshots
        Xvfb :1 -screen 0 1024x768x16 &> xvfb.log &
        export DISPLAY=:1

        # Run the IGV screenshot script with the provided inputs
        python3 /opt/IGV_Linux_2.18.2/make_igv_screenshot.py \
          ~{aligned_bam_hap1} ~{aligned_bam_hap2} ~{alignments} \
          ~{aligned_bam_hap1_bai} ~{aligned_bam_hap2_bai} ~{alignments_bai} \
          -r ~{bed_file} \
          -ht ~{image_height} \
          -bin /opt/IGV_Linux_2.18.2/igv.sh \
          -mem ~{memory_mb} \
          --fasta_file ~{fasta_file} \
          --fasta_fai ~{fasta_fai} \
          --sample_name ~{sample_name}
    }

    runtime {
        docker: docker_image
        memory: "~{memory_mb} MB"
        cpu: 2
        disks: "local-disk ~{disk_gb} HDD"
    }

    output {
        Array[File] snapshots = glob("IGV_Snapshots/*.png")
    }
}
