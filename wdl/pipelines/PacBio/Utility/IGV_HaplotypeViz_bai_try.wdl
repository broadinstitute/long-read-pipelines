version 1.0

workflow IGV_HaplotypeViz {
  input {
    File aligned_bam_hap1
    File aligned_bam_hap2
    File alignments_bam
    File ref_fasta
    File targeted_bed_file
    String sample_name
    String disk_type = "SSD"
    String output_dir
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
      output_dir = output_dir
  }

  output {
    Array[File] final_screenshots = FinalizeToGCS.uploaded_files
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
    mkdir -p /cromwell_root/output/IGV_Snapshots && chmod 777 /cromwell_root/output/IGV_Snapshots

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
    Array[File] screenshots = glob("/cromwell_root/output/IGV_Snapshots/*.png")
  }

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
    set -euxo pipefail

    # Ensure the output directory exists and is properly formatted
    gcs_output_dir=$(echo ~{output_dir} | sed 's:/*$::')

    # Copy all screenshots to Google Cloud Storage
    for file in ~{sep=' ' screenshots}; do
      gsutil cp $file $gcs_output_dir/
    done
  >>>

  output {
    Array[File] uploaded_files = glob("~{output_dir}/*.png")
  }

  runtime {
    cpu:        1
    memory:     "4 GiB"
    disks:      "local-disk 10 SSD"
    preemptible: 1
    maxRetries: 2
    docker:     "google/cloud-sdk:slim"
  }
}
