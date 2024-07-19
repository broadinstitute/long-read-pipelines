version 1.0

import "../../structs/Structs.wdl"

task processWithTRGT {
  input {
    File input_bam
    File input_bam_bai
    String basename
    File ref_fasta
    File ref_fasta_index
    String repeatCatalog
    Int is_female # 1=female, 0=male
    Int cpuCores = 16

    RuntimeAttr? runtime_attr_override

  }
  
  meta {
    description: "Uses TRGT to size TRs in a bam file."
  }

  command <<<
    set -euo pipefail
    karyotype="XY"
    if [[ ~{is_female} -eq 1 ]]; then karyotype="XX"; fi 
    trgt genotype --genome ~{ref_fasta} --repeats ~{repeatCatalog} --reads ~{input_bam} --threads ~{cpuCores} --output-prefix ~{basename}_trgt --karyotype ${karyotype}
    bcftools sort -Ob -o ~{basename}_trgt_sorted.vcf.gz ~{basename}_trgt.vcf.gz
    mv ~{basename}_trgt_sorted.vcf.gz ~{basename}_trgt.vcf.gz
    bcftools index -t ~{basename}_trgt.vcf.gz

  >>>

  output {
    File trgt_output_vcf = "~{basename}_trgt.vcf.gz"
    File trgt_output_vcf_idx = "~{basename}_trgt.vcf.gz.tbi"
    File trgt_output_bam = "~{basename}_trgt.spanning.bam"
  }
  
  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          cpuCores,
      mem_gb:             16,
      disk_gb:            500,
      boot_disk_gb:       10,
      preemptible_tries:  3,
      max_retries:        1,
      docker:             "us.gcr.io/broad-dsp-lrma/lr-trgt:1.0.0"
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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
