version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/VariantCalling/ONTPepper.wdl" as ONTPepper
import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/QC/SampleLevelAlignedMetrics.wdl" as COV

workflow ONTWholeGenomePlasmid {

    meta {
        description: "A workflow that performs single sample variant calling on Oxford Nanopore reads from one or more flow cells. The workflow merges multiple flowcells into a single BAM prior to variant calling."
    }
    parameter_meta {
        aligned_bam:       "GCS path to aligned BAM files"
        aligned_bai:       "GCS path to aligned BAM file indices"
        participant_name:  "name of the participant from whom these samples were obtained"

        ref_fasta:         "FASTA"
        ref_fasta_fai:     "FASTA fai"
        gcs_out_root_dir:  "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        File aligned_bam
        File aligned_bai
        String participant_name

        File ref_fasta
        File ref_fasta_fai
        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTWholeGenomePlasmid/~{participant_name}"

    String dir = outdir + "/alignments"

    call ONTPepper.Pepper {
        input:
            bam           = aligned_bam,
            bai           = aligned_bai,
            ref_fasta     = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            threads       = 8,
            memory        = 64
    }

    output {
        File dvp_vcf = Pepper.VCF
        File dvp_tbi = Pepper.VCF_tbi
    }
}
