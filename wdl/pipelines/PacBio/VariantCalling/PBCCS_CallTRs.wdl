version 1.0

import "../../../tasks/VariantCalling/TRGT.wdl" as TRGT
import "../../../structs/Structs.wdl"

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
    String repeatCatalog = "adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0.bed"

    #Optional runtime arguments
    RuntimeAttr? runtime_attr_override
  }

  call TRGT.processWithTRGT as processWithTRGT {
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
