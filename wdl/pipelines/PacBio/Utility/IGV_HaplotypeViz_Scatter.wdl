version 1.0

import "../../../structs/Structs.wdl"

workflow IGV_HaplotypeViz_Scatter {
  input {
    # BED files containing regions to screenshot; 4th column can optionally be SVID
    Array[File] beds
    Array[String]? run_names

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

    # Docker image for IGV headless tasks
    String igv_docker
  }

  scatter (i in range(length(beds))) {
    String sample_combined = sample_id + "_combined"

    # Run IGV for both BAM and FASTA visualization for both haplotypes (Hap1 and Hap2)
    call RunIGVHeadlessCombined {
      input:
        bam_hap1=bam_hap1,
        bai_hap1=bai_hap1,
        bam_hap2=bam_hap2,
        bai_hap2=bai_hap2,
        fasta_hap1=haplotig_fasta_hap1,
        fasta_hap2=haplotig_fasta_hap2,
        bed=beds[i],
        sample_id=sample_combined,
        ref_fasta=ref_fasta,
        ref_fai=ref_fai,
        igv_docker=igv_docker
    }
  }

  output {
    Array[File] igv_screenshots_combined = RunIGVHeadlessCombined.igv_screenshot
  }
}

task RunIGVHeadlessCombined {
  input {
    File bam_hap1         # BAM file for Haplotype 1
    File bai_hap1         # BAI file for Haplotype 1
    File bam_hap2         # BAM file for Haplotype 2
    File bai_hap2         # BAI file for Haplotype 2
    File fasta_hap1       # FASTA file for Haplotype 1
    File fasta_hap2       # FASTA file for Haplotype 2
    File bed              # BED file containing regions to visualize (3 or 4 columns allowed)
    String sample_id      # Sample ID for naming outputs
    File ref_fasta        # Reference genome used for alignment
    File ref_fai          # Index for the reference genome
    String igv_docker     # Docker image for running IGV headless
    Int? records_per_shard # Optional: Parallelization parameter for large datasets
  }

  command <<<
    # Running IGV headless mode to take screenshots for both BAM and FASTA files for both haplotypes
    igv.sh \
    -b ~{bam_hap1},~{bam_hap2} \
    -i ~{bai_hap1},~{bai_hap2} \
    -g ~{fasta_hap1},~{fasta_hap2} \
    -bed ~{bed} \
    -name hap1_bam,hap2_bam,hap1_fasta,hap2_fasta \
    -o ~{sample_id}.igv_screenshot.png
  >>>

  output {
    File igv_screenshot = "~{sample_id}.igv_screenshot.png"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/igv_docker:v952024"  # Updated IGV docker image
    memory: "8G"
    cpu: "2"
    disks: "local-disk 10 HDD"
  }
}
