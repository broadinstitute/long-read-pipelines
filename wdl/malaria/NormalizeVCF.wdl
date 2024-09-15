version 1.0

workflow NormalizeVCF {
  # Define inputs in an input block
  input {
    File input_vcf
    File reference_fa
  }

  # Step 1: Remove HAPCOMP field
  call RemoveHAPCOMP {
    input:
      input_vcf = input_vcf
  }

  # Step 2: Remove HAPDOM field
  call RemoveHAPDOM {
    input:
      input_vcf = RemoveHAPCOMP.output_vcf
  }

  # Step 3: Normalize VCF
  call NormalizeVCFFile {
    input:
      input_vcf = RemoveHAPDOM.output_vcf,
      reference_fa = reference_fa
  }

  # Output
  output {
    File normalized_vcf = NormalizeVCFFile.output_vcf
  }
}

task RemoveHAPCOMP {
  input {
    File input_vcf
  }

  command {
    bcftools annotate -x 'INFO/HAPCOMP' ~{input_vcf} | bgzip -c > output.no_hapcomp.vcf.gz
  }

  output {
    File output_vcf = "output.no_hapcomp.vcf.gz"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/bcftools_htslib:v9152024"
    memory: "4G"
    cpu: 1
  }
}

task RemoveHAPDOM {
  input {
    File input_vcf
  }

  command {
    bcftools annotate -x 'INFO/HAPDOM' ~{input_vcf} | bgzip -c > output.no_hapdom.vcf.gz
  }

  output {
    File output_vcf = "output.no_hapdom.vcf.gz"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/bcftools_htslib:v9152024"
    memory: "4G"
    cpu: 1
  }
}

task NormalizeVCFFile {
  input {
    File input_vcf
    File reference_fa
  }

  command {
    bcftools norm -m -any --atom-overlaps . -f ~{reference_fa} ~{input_vcf} | bgzip -c > output.norm.vcf.gz
  }

  output {
    File output_vcf = "output.norm.vcf.gz"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/bcftools_htslib:v9152024"
    memory: "4G"
    cpu: 1
  }
}
