version 1.0

import "../../../tasks/Utility/AnnotationUtils.wdl" as AU

workflow AnnotateVariantContext {

    meta {
        description: "A workflow that annotates a VCF with sequence context information"
    }

    parameter_meta {
        vcf_gz:       "The VCF to annotate"
        vcf_gz_tbi:   "The VCF tabix index"
        ref_fa:       "Reference FASTA file"
    }

    input {
        File vcf_gz
        File vcf_gz_tbi
        File ref_fasta
    }

    call AU.AnnotateVCF {
        input:
            vcf_gz = vcf_gz,
            vcf_gz_tbi = vcf_gz_tbi,
            ref_fasta = ref_fasta
    }

    output {
        File annotated_vcf_gz = AnnotateVCF.annotated_vcf_gz 
        File annotated_vcf_gz_tbi = AnnotateVCF.annotated_vcf_gz_tbi
    }
}
