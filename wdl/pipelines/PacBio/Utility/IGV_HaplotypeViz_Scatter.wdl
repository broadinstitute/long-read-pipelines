version 1.0

import "../../../structs/Structs.wdl"

workflow IGV_HaplotypeViz_Scatter {
  input {
    # BED files containing regions to screenshot; 4th column can optionally be SVID
    Array[File] beds
    Array[String] run_names

    # BAM and BAI files from align_asm workflow for alignment visualization
    File bam_hap1
    File bai_hap1
    File bam_hap2
    File bai_hap2

    # FASTA files from PBAssembleWithHifiasm or bam_to_contig workflow for sequence visualization
    File haplotig_fasta_hap1
    File haplotig_fasta_hap2

    # Reference corresponding to read alignments for BAM files
    File ref_fasta
    File ref_fai

    # Sample id and prefix for output filenames
    String sample_id

    # Number of records per shard for parallelization
    Int? records_per_shard

    # Docker images for Linux and IGV headless tasks
    String linux_docker
    String igv_docker
  }

  scatter (i in range(length(beds))) {
    String sample_w_hap1 = sample_id + "_hap1"
    String sample_w_hap2 = sample_id + "_hap2"

    # Run IGV for BAM alignments Haplotype 1 (H1)
    call RunIGVHeadless as IGV_Hap1 {
      input:
        bam_or_cram=bam_hap1,
        bam_or_cram_index=bai_hap1,
        bed=beds[i],
        sample_id=sample_w_hap1,
        ref_fasta=ref_fasta,
        ref_fai=ref_fai,
        igv_docker=igv_docker
    }

    # Run IGV for BAM alignments Haplotype 2 (H2)
    call RunIGVHeadless as IGV_Hap2 {
      input:
        bam_or_cram=bam_hap2,
        bam_or_cram_index=bai_hap2,
        bed=beds[i],
        sample_id=sample_w_hap2,
        ref_fasta=ref_fasta,
        ref_fai=ref_fai,
        igv_docker=igv_docker
    }

    # For sequence visualization, use FASTA for haplotigs Haplotype 1 (H1)
    call RunIGVHeadless as IGV_Seq_Hap1 {
      input:
        bam_or_cram=haplotig_fasta_hap1,
        bam_or_cram_index=bai_hap1,  # Index may not be necessary for FASTA
        bed=beds[i],
        sample_id=sample_w_hap1,
        ref_fasta=ref_fasta,  # Reference may not be needed for FASTA visualization
        ref_fai=ref_fai,
        igv_docker=igv_docker
    }

    # For sequence visualization, use FASTA for haplotigs Haplotype 2 (H2)
    call RunIGVHeadless as IGV_Seq_Hap2 {
      input:
        bam_or_cram=haplotig_fasta_hap2,
        bam_or_cram_index=bai_hap2,
        bed=beds[i],
        sample_id=sample_w_hap2,
        ref_fasta=ref_fasta,
        ref_fai=ref_fai,
        igv_docker=igv_docker
    }
  }

  output {
    Array[File] igv_screenshots_hap1 = IGV_Hap1.igv_screenshot
    Array[File] igv_screenshots_hap2 = IGV_Hap2.igv_screenshot
    Array[File] igv_screenshots_seq_hap1 = IGV_Seq_Hap1.igv_screenshot
    Array[File] igv_screenshots_seq_hap2 = IGV_Seq_Hap2.igv_screenshot
  }
}

task RunIGVHeadless {
  input {
    File bam_or_cram       # BAM/CRAM or FASTA file for visualization
    File bam_or_cram_index # Index file for BAM/CRAM
    File bed               # BED file containing regions to visualize (3 or 4 columns allowed)
    String sample_id       # Sample ID for naming outputs
    File ref_fasta         # Reference genome used for alignment
    File ref_fai           # Index for the reference genome
    String igv_docker      # Docker image for running IGV headless
    Int? records_per_shard # Optional: Parallelization parameter for large datasets
  }

  command <<<
    # Running IGV headless mode to take screenshots for each region in the BED file
    igv.sh \
    -b ~{bam_or_cram} \
    -i ~{bam_or_cram_index} \
    -g ~{ref_fasta} \
    -bed ~{bed} \
    -o ~{sample_id}.igv_screenshot.png
  >>>

  output {
    File igv_screenshot = "~{sample_id}.igv_screenshot.png"
  }

  runtime {
    docker: "~{igv_docker}"
    memory: "8G"
    cpu: "2"
    disks: "local-disk 10 HDD"
  }
}
