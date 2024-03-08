version 1.0

import "../../../tasks/Utility/AnnotationUtils.wdl" as AU
import "../../../tasks/Utility/VariantUtils.wdl" as VU

workflow AnnotateVariantContext {

    meta {
        description: "A workflow that annotates a VCF with sequence context information"
    }

    parameter_meta {
        vcf_gz:    "The VCF to annotate"
        vcf_tbi:   "The VCF tabix index"
        ref_fasta: "Reference FASTA file"
    }

    input {
        File vcf_gz
        File vcf_tbi
        File ref_fasta
    }

    call AU.AnnotateVCF {
        input:
            vcf_gz = vcf_gz,
            vcf_tbi = vcf_tbi,
            ref_fasta = ref_fasta
    }

    call VU.ZipAndIndexVCF { input: vcf = AnnotateVCF.annotated_vcf }

    output {
        File annotated_vcf_gz = ZipAndIndexVCF.vcfgz
        File annotated_tbi = ZipAndIndexVCF.tbi
    }
}
