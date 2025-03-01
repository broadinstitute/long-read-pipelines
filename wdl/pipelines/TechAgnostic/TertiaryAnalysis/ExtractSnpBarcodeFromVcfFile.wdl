version 1.0

import "../../../structs/Structs.wdl"
import "../../../tasks/Utility/VariantUtils.wdl" as VariantUtils

workflow ExtractSnpBarcodeFromVcfFile  {

    meta {
        author: "Jonn Smith"
        description: "A workflow that extracts a barcode from a VCF file using a haplotype database file."
    }

    parameter_meta {
        vcf: "The VCF file to extract the barcode from."
        vcf_index: "The index file for the VCF file."
        haplotype_database_file: "The haplotype database file."
        ref_fasta: "The reference fasta file."
        ref_fasta_fai: "The index file for the reference fasta file."
        ref_dict: "The dictionary file for the reference fasta file."
        prefix: "The prefix for the output files."
        runtime_attr_override: "The runtime attributes to override."
    } 

    input {
        File vcf
        File vcf_index

        File haplotype_database_file

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String prefix = "barcode"

        RuntimeAttr? runtime_attr_override
    }
    
    call VariantUtils.ExtractFingerprintAndBarcode as t01_ExtractFingerprintAndBarcode {
        input:
            vcf = vcf,
            vcf_index = vcf_index,
            haplotype_database_file = haplotype_database_file,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_dict = ref_dict,
            prefix = prefix,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File output_vcf = t01_ExtractFingerprintAndBarcode.output_vcf
        String barcode = t01_ExtractFingerprintAndBarcode.barcode
        File barcode_file = t01_ExtractFingerprintAndBarcode.barcode_file
    }
}
