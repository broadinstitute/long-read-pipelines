version 1.0

workflow runTRGT {

  meta {
    description: "Uses TRGT to size TRs in a bam file."
  }

  input {
    File input_bam
    File input_bam_bai
    String basename = basename(input_bam, ".bam")
    String output_gs_path
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File repeatCatalog = "https://zuchnerlab.s3.amazonaws.com/RepeatExpansions/TRGT/adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0.bed"

    #Optional runtime arguments
    RuntimeAttr? runtime_attr_override
  }

  call processWithTRGT {
    input:
      input_bam = input_bam,
      input_bam_bai = input_bam_bai,
      basename = basename,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      repeatCatalog = repeatCatalog,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File trgt_output_vcf = processWithTRGT.trgt_output_vcf
    File trgt_output_bam = processWithTRGT.trgt_output_bam
  }
}

task processWithTRGT {
  input {
    File input_bam
    File input_bam_bai
    String basename
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File repeatCatalog

    RuntimeAttr? runtime_attr_override

  }
  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          4,
      mem_gb:             16,
      disk_gb:            500,
      boot_disk_gb:       10,
      preemptible_tries:  3,
      max_retries:        1,
      docker:             "public.ecr.aws/s5z5a3q9/lr-trgt:0.4.0"
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  meta {
    description: "Uses TRGT to size TRs in a bam file."
  }

  command <<<
    set -euo pipefail
    trgt --genome ~{ref_fasta} --repeats ~{repeatCatalog} --reads ~{input_bam} --threads ~{runtime_attr.cpu_cores} --output-prefix ~{basename}_trgt

  >>>

  output {
    File trgt_output_vcf = "~{basename}_trgt.vcf.gz"
    File trgt_output_bam = "~{basename}_trgt.spanning.bam"
  }

  runtime {
      cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
      memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
      disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
      bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
      preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
      maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
      docker:                 select_first([runtime_attr.docker,            default_attr.docker])
  }
}