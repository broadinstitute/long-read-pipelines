version 1.0

workflow igv_screenshot_workflow {
  input {
    File aligned_bam_hap1          # BAM file for haplotype 1
    File aligned_bam_hap2          # BAM file for haplotype 2
    File alignments_bam            # Total alignments BAM file
    File ref_fasta                 # Reference FASTA file
    File targeted_bed_file         # BED file with regions of interest
    String sample_name             # Name for the sample (used in output naming)
    String disk_type = "SSD"       # Default disk type
    String gcs_output_dir          # GCS directory to copy outputs
  }

  call GenerateIgvScreenshots {
    input:
      aligned_bam_hap1 = aligned_bam_hap1,
      aligned_bam_hap2 = aligned_bam_hap2,
      alignments_bam = alignments_bam,
      ref_fasta = ref_fasta,
      targeted_bed_file = targeted_bed_file,
      sample_name = sample_name,
      disk_type = disk_type
  }

  call FinalizeToGCS {
    input:
      screenshots = GenerateIgvScreenshots.screenshots,
      output_dir = gcs_output_dir
  }

  output {
    Array[File] screenshot_files = GenerateIgvScreenshots.screenshots
  }
}

task GenerateIgvScreenshots {
  input {
    File aligned_bam_hap1
    File aligned_bam_hap2
    File alignments_bam
    File ref_fasta
    File targeted_bed_file
    String sample_name
    String disk_type
  }

  command <<<
    # Ensure the snapshots directory exists and set permissions
    mkdir -p /output/IGV_Snapshots && chmod 777 /output/IGV_Snapshots

    # Start a virtual frame buffer to allow IGV to render
    Xvfb :1 -screen 0 1024x768x16 & export DISPLAY=:1

    # Run the Python script to generate IGV snapshots
    python3 /opt/IGV_Linux_2.18.2/make_igv_screenshot.py \
      ~{aligned_bam_hap1} ~{aligned_bam_hap2} ~{alignments_bam} \
      --fasta_file ~{ref_fasta} \
      --sample_name ~{sample_name} \
      -r ~{targeted_bed_file} \
      -ht 500 \
      -bin /opt/IGV_Linux_2.18.2/igv.sh \
      -mem 4000

  >>>

  output {
    Array[File] screenshots = glob("/output/IGV_Snapshots/*.png")
  }

  # Calculate dynamic disk size based on the size of BAM files
  Int disk_size = ceil(size(aligned_bam_hap1, "GiB")) + ceil(size(aligned_bam_hap2, "GiB")) + ceil(size(alignments_bam, "GiB")) + 10

  runtime {
    cpu:        4
    memory:     "8 GiB"
    disks:      "local-disk " + disk_size + " " + disk_type
    preemptible: 2
    maxRetries: 1
    docker:     "us.gcr.io/broad-dsp-lrma/igv_screenshot_docker:v982024"
  }
}

task FinalizeToGCS {
  input {
    Array[File] screenshots
    String output_dir
  }

  command <<<
    # Copy the output PNG files to Google Cloud Storage
    for file in ~{sep=' ' screenshots}; do
      gsutil cp "$file" "~{output_dir}/"
    done
  >>>

  output {
    Array[File] uploaded_files = screenshots
  }

  runtime {
    cpu:        1
    memory:     "2 GiB"
    disks:      "local-disk 10 HDD"
    preemptible: 2
    maxRetries: 1
    docker:     "gcr.io/google-containers/toolbox:latest"
  }
}
