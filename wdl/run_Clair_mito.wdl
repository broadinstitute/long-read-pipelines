version 1.0

import "tasks/Clair_mito.wdl" as Clair_mito
import "tasks/Structs.wdl"

workflow run_Clair_mito {
    input {
        File bam
        File bai
        File ref_fasta
        File ref_fasta_fai
        File? sites_vcf
        File? sites_vcf_tbi
        String? chr
        String preset
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta{
        bam:        "GCS path to raw subread bam"
        bai:        "index for bam file"
        locus:      "genomic locus to select"
        prefix:     "prefix for output bam and bai file names"
        ref_fasta:  "chrM reference fasta"
        map_preset: "preset to be used for minimap2 parameter '-x'"
        ref_fai:    "index of fa"

    }

    call Clair_mito.Clair as Clair_mito_call {
        input:
            bam = bam,
            bai = bai,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            preset = preset
     }

    output {
            File pileup_vcf = Clair_mito_call.pileup_vcf
            File pileup_vcf_tbi = Clair_mito_call.pileup_vcf_tbi
            File full_alignment_vcf = Clair_mito_call.full_alignment_vcf
            File full_alignment_tbi = Clair_mito_call.full_alignment_tbi

            File vcf = Clair_mito_call.vcf
            File vcf_tbi = Clair_mito_call.vcf_tbi
            File gvcf = Clair_mito_call.gvcf
            File gvcf_tbi = Clair_mito_call.gvcf_tbi
        }

}