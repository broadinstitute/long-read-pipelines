version 1.0

workflow NormalizeVCF {
  input {
    File input_vcf
    File reference_fa
    String sample_name  # New input for sample name
    Int disk_size_gb = 20
    Int memory_gb = 4
    Int cpu_cores = 1
  }

  call RemoveHAPCOMP {
    input:
      input_vcf = input_vcf,
      disk_size_gb = disk_size_gb,
      memory_gb = memory_gb,
      cpu_cores = cpu_cores
  }

  call RemoveHAPDOM {
    input:
      input_vcf = RemoveHAPCOMP.output_vcf,
      disk_size_gb = disk_size_gb,
      memory_gb = memory_gb,
      cpu_cores = cpu_cores
  }

  call NormalizeVCFFile {
    input:
      input_vcf = RemoveHAPDOM.output_vcf,
      reference_fa = reference_fa,
      sample_name = sample_name,
      disk_size_gb = disk_size_gb,
      memory_gb = memory_gb,
      cpu_cores = cpu_cores
  }

  output {
    File normalized_vcf = NormalizeVCFFile.output_vcf
  }
}

task RemoveHAPCOMP {
  input {
    File input_vcf
    Int disk_size_gb
    Int memory_gb
    Int cpu_cores
  }

  command {
    bcftools annotate -x 'INFO/HAPCOMP' ~{input_vcf} | bgzip -c > output.no_hapcomp.vcf.gz
  }

  output {
    File output_vcf = "output.no_hapcomp.vcf.gz"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/bcftools_htslib:v9152024"
    memory: "~{memory_gb}G"
    cpu: "~{cpu_cores}"
    disks: "local-disk ~{disk_size_gb} HDD"
  }
}

task RemoveHAPDOM {
  input {
    File input_vcf
    Int disk_size_gb
    Int memory_gb
    Int cpu_cores
  }

  command {
    bcftools annotate -x 'INFO/HAPDOM' ~{input_vcf} | bgzip -c > output.no_hapdom.vcf.gz
  }

  output {
    File output_vcf = "output.no_hapdom.vcf.gz"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/bcftools_htslib:v9152024"
    memory: "~{memory_gb}G"
    cpu: "~{cpu_cores}"
    disks: "local-disk ~{disk_size_gb} HDD"
  }
}

task NormalizeVCFFile {
  input {
    File input_vcf
    File reference_fa
    String sample_name  # Input for the sample name
    Int disk_size_gb
    Int memory_gb
    Int cpu_cores
  }

  command {
    bcftools norm -m -any --atom-overlaps . -f ~{reference_fa} ~{input_vcf} | bgzip -c > ~{sample_name}.norm.vcf.gz
  }

  output {
    File output_vcf = "~{sample_name}.norm.vcf.gz"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/bcftools_htslib:v9152024"
    memory: "~{memory_gb}G"
    cpu: "~{cpu_cores}"
    disks: "local-disk ~{disk_size_gb} SSD"
  }
}
