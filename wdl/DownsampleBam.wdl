version 1.0

workflow DownsampleBam {

  meta {
    description: "Uses Picard DownsampleSam to downsample a Bam file."
  }

  input {
    File input_bam
    File input_bam_bai
    String basename = basename(input_bam, ".bam")
    String output_gs_path
    Int desiredCoverage
    Float currentCoverage

    #Optional runtime arguments
    Int? preemptible_tries

  }

  parameter_meta {
    basename: "Input is a string specifying the sample name which will be used to locate the file on gs."
    desiredCoverage: "Input is an integer of the desired approximate coverage in the output bam file."
    output_gs_path: "GCS folder path name to which to save the output files."
    downsampled_bam: "Output is a bam file downsampled to the specified mean coverage."
    downsampled_bai: "Output is the index file for a bam file downsampled to the specified mean coverage."
  }

  call downsampleBam {
    input:
      input_bam = input_bam,
      input_bam_bai = input_bam_bai,
      basename = basename,
      output_gs_path = output_gs_path,
      desiredCoverage = desiredCoverage,
      currentCoverage = currentCoverage,
      preemptible_tries = preemptible_tries
  }

  output {
    String downsampled_bam = downsampleBam.downsampled_bam
    String downsampled_bai = downsampleBam.downsampled_bai
  }
}

task downsampleBam {
  input {
    File input_bam
    File input_bam_bai
    String basename
    String output_gs_path
    Int desiredCoverage
    Float currentCoverage
    Float scalingFactor = desiredCoverage / currentCoverage

    Int? preemptible_tries
  }

  meta {
    description: "Uses Picard to downsample to desired coverage based on provided estimate of coverage."
  }
  parameter_meta {
    basename: "Input is a string specifying the sample name which will be used to locate the file on gs."
    output_gs_path: "GCS folder path name to which to save the output files."
    downsampled_bam: "Output is a bam file downsampled to the specified mean coverage."
    downsampled_bai: "Output is the index file for a bam file downsampled to the specified mean coverage."
    desiredCoverage: "Input is an integer of the desired approximate coverage in the output bam file."
  }
  command <<<
    set -eo pipefail
    gatk DownsampleSam -I ~{input_bam} -O ~{basename}_~{desiredCoverage}x.bam -R 7 -P ~{scalingFactor} -S ConstantMemory --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true
    gsutil cp "~{basename}_~{desiredCoverage}x.bam" "~{output_gs_path}"
    gsutil cp "~{basename}_~{desiredCoverage}x.bai" "~{output_gs_path}"


  >>>
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "8 GB"
    cpu: "2"
    disks: "local-disk 500 HDD"
    docker: "us.gcr.io/broad-gatk/gatk"
  }
  output {
    String downsampled_bam = "~{output_gs_path}~{basename}_~{desiredCoverage}x.bam"
    String downsampled_bai = "~{output_gs_path}~{basename}_~{desiredCoverage}x.bai"
  }
}