version 1.0

import "ONTPepper.wdl"

import "../../tasks/Utility/VariantUtils.wdl"
import "../../tasks/Utility/Utils.wdl"
import "../../tasks/Alignment/WhatsHap.wdl"

workflow Run {
    meta {
        desciption:
        "Runs Clair3 on the input (sharded) BAM."
    }
    parameter_meta {
        how_to_shard_wg_for_calling: "An array of the BAM's shard; each element is assumed to be a tuple of (ID for the shard, (BAM of the shard, BAI of the shard))"
        prefix: "Prefix for output files"
        model_for_pepper_margin_dv: "refer to https://github.com/kishwarshafin/pepper for appropriate values"
    }

    input {
        Array[Pair[String, Pair[File, File]]] how_to_shard_wg_for_calling
        String prefix
        String model_for_pepper_margin_dv

        Map[String, String] ref_map

        # optimization
        Int dv_threads
        Int dv_memory
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    output {
        File legacy_ont_dvp_g_vcf = MergePEPPERGVCFs.vcf
        File legacy_ont_dvp_g_tbi = MergePEPPERGVCFs.tbi
        File legacy_ont_dvp_phased_vcf = MergePEPPERPhasedVCFs.vcf
        File legacy_ont_dvp_phased_tbi = MergePEPPERPhasedVCFs.tbi
        File legacy_ont_dvp_haplotagged_bam = MergePEPPERHapTaggedBam.merged_bam
        File legacy_ont_dvp_haplotagged_bai = MergePEPPERHapTaggedBam.merged_bai
        File legacy_ont_dvp_phased_vcf_stats_tsv = ONTPhaseStatsLegacy.stats_tsv
        File legacy_ont_dvp_phased_vcf_stats_gtf = ONTPhaseStatsLegacy.stats_gtf
    }

    scatter (triplet in how_to_shard_wg_for_calling) {
        call ONTPepper.Pepper {
            input:
                bam           = triplet.right.left,
                bai           = triplet.right.right,
                ref_fasta     = ref_map['fasta'],
                ref_fasta_fai = ref_map['fai'],
                model         = model_for_pepper_margin_dv,
                threads       = select_first([dv_threads]),
                memory        = select_first([dv_memory]),
                zones         = zones
        }
    }

    String pepper_prefix = prefix + ".PEPPER-Margin-DeepVariant"

    call VariantUtils.MergeAndSortVCFs as MergePEPPERGVCFs {
        input:
            vcfs     = Pepper.gVCF,
            prefix   = pepper_prefix + ".g",
            ref_fasta_fai = ref_map['fai']
    }

    call VariantUtils.MergeAndSortVCFs as MergePEPPERPhasedVCFs {
        input:
            vcfs     = Pepper.phasedVCF,
            prefix   = pepper_prefix + ".phased",
            ref_fasta_fai = ref_map['fai']
    }

    call Utils.MergeBams as MergePEPPERHapTaggedBam {
        input:
            bams     = Pepper.hap_tagged_bam,
            prefix   = prefix +  ".MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged"
    }

    call WhatsHap.Stats as ONTPhaseStatsLegacy {
        input:
        phased_vcf=MergePEPPERPhasedVCFs.vcf, phased_tbi=MergePEPPERPhasedVCFs.tbi
    }
}
