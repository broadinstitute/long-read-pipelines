version 1.0

import "Utils.wdl"

task HaplotypeCaller_GATK4_VCF {
  input {
    File input_bam
    File input_bam_index
    File interval_list
    String vcf_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? contamination
    Boolean make_gvcf
    Boolean make_bamout
    Int preemptible_tries
    Int hc_scatter
    Boolean run_dragen_mode_variant_calling = false
    Boolean use_dragen_hard_filtering = false
    Boolean use_spanning_event_genotyping = true
    File? dragstr_model
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
    Int memory_multiplier = 1
  }

  Int memory_size_mb = ceil(8000 * memory_multiplier)

  String output_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
  String output_file_name = vcf_basename + output_suffix

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(((size(input_bam, "GiB") + 30) / hc_scatter) + ref_size) + 20

  String bamout_arg = if make_bamout then "-bamout ~{vcf_basename}.bamout.bam" else ""

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command <<<
    set -e
    # We need at least 1 GB of available memory outside of the Java heap in order to execute native code, thus, limit
    # Java's memory by the total memory minus 1 GB. We need to compute the total memory as it might differ from
    # memory_size_gb because of Cromwell's retry with more memory feature.
    # Note: In the future this should be done using Cromwell's ${MEM_SIZE} and ${MEM_UNIT} environment variables,
    #       which do not rely on the output format of the `free` command.
    available_memory_mb=$(free -m | awk '/^Mem/ {print $2}')
    let java_memory_size_mb=available_memory_mb-1024
    echo Total available memory: ${available_memory_mb} MB >&2
    echo Memory reserved for Java: ${java_memory_size_mb} MB >&2

    gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
      HaplotypeCaller \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -L ~{interval_list} \
      -O ~{output_file_name} \
      -contamination ~{default=0 contamination} \
      -G StandardAnnotation -G StandardHCAnnotation ~{true="-G AS_StandardAnnotation" false="" make_gvcf} \
      ~{true="--dragen-mode" false="" run_dragen_mode_variant_calling} \
      ~{false="--disable-spanning-event-genotyping" true="" use_spanning_event_genotyping} \
      ~{if defined(dragstr_model) then "--dragstr-params-path " + dragstr_model else ""} \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      ~{true="-ERC GVCF" false="" make_gvcf} \
      ~{bamout_arg}

    # Cromwell doesn't like optional task outputs, so we have to touch this file.
    touch ~{vcf_basename}.bamout.bam
  >>>

  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "~{memory_size_mb} MiB"
    cpu: "2"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_vcf = "~{output_file_name}"
    File output_vcf_index = "~{output_file_name}.tbi"
    File bamout = "~{vcf_basename}.bamout.bam"
  }
}