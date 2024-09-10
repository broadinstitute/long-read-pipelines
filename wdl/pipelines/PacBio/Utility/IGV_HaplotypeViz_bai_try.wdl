version 1.0

workflow igv_screenshot_workflow {
  input {
    File aligned_bam_hap1          # BAM file for haplotype 1
    File aligned_bam_hap2          # BAM file for haplotype 2
    File alignments_bam            # Total alignments BAM file
    File ref_fasta                 # Reference FASTA file
    File targeted_bed_file         # BED file with regions of interest
    String sample_name             # Sample name for naming convention
    Int image_height = 500         # Height of IGV track, default to 500
    Int memory_mb = 4000           # Memory for IGV, default to 4000MB
  }

  call make_igv_screenshot {
    input:
      aligned_bam_hap1 = aligned_bam_hap1,
      aligned_bam_hap2 = aligned_bam_hap2,
      alignments_bam = alignments_bam,
      ref_fasta = ref_fasta,
      targeted_bed_file = targeted_bed_file,
      sample_name = sample_name,
      image_height = image_height,
      memory_mb = memory_mb
  }

  output {
    Array[File] pngs = make_igv_screenshot.pngs  # Collect all generated PNG files
  }
}

task make_igv_screenshot {
  input {
    File aligned_bam_hap1          # BAM file for haplotype 1
    File aligned_bam_hap2          # BAM file for haplotype 2
    File alignments_bam            # Total alignments BAM file
    File ref_fasta                 # Reference FASTA file
    File targeted_bed_file         # BED file with regions of interest
    String sample_name             # Sample name for naming convention
    Int image_height               # Height of IGV track
    Int memory_mb                  # Memory for IGV
  }

  command {
    # Create output directory for snapshots
    mkdir -p IGV_Snapshots
    
    # Start a virtual framebuffer (Xvfb) to allow IGV to render without display
    Xvfb :1 -screen 0 1024x768x16 &> xvfb.log &
    export DISPLAY=:1

    # Run the Python script to generate IGV screenshots
    python3 /opt/IGV_Linux_2.18.2/make_igv_screenshot.py \
      ~{aligned_bam_hap1} ~{aligned_bam_hap2} ~{alignments_bam} \
      -r ~{targeted_bed_file} \
      -ht ~{image_height} \
      -bin /opt/IGV_Linux_2.18.2/igv.sh \
      -mem ~{memory_mb} \
      --fasta_file ~{ref_fasta} \
      --sample_name ~{sample_name}
  }

  output {
    # Capture all generated PNG snapshot files
    Array[File] pngs = glob("IGV_Snapshots/*.png")
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/igv_screenshot_docker:v982024"
    memory: "~{memory_mb} MB"
    cpu: 2
    disks: "local-disk 50 HDD"  # Specify disk size if needed
  }
}
