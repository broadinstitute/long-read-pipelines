version 1.0

import "../Utility/Utils.wdl"
import "../Utility/VariantUtils.wdl"

import "ShardWholeGenome.wdl"

import "DeepVariant.wdl"
import "MarginPhase.wdl" as Margin
import "../Alignment/WhatsHap.wdl"

import "Clair.wdl" as Clair3
import "../../deprecated/tasks/ONTPepper.wdl"

workflow Work {
    meta {
        description: "Call small variants using reads-based methods (i.e. not for assembly-contig-based methods)."
    }
    parameter_meta {
        prefix: "Prefix for output files"
        per_chr_bam_bai_and_id: "WGS bam sharded per chromosome/contig."
        is_ont: "If the input data is ONT"
        is_r10_4_pore_or_later: "If the ONT input data is generated on R10.4 simples/duplex pores."
        model_for_dv_andor_pepper: "Model string to be used on DV or the PEPPER-Margin-DeepVariant toolchain. Please refer to their github pages for accepted values."
        ref_scatter_interval_list_locator: "A file holding paths to interval_list files, used for custom sharding the of the input BAM; when not provided, will shard WG by contig (possibly slower)"
        ref_scatter_interval_list_ids: "A file that gives short IDs to the interval_list files; when not provided, will shard WG by contig (possibly slower)"
        use_gpu: "Use GPU acceleration for DV (or PEPPER) or not"
    }
    input {
        # sample info
        File bam
        File bai
        String prefix

        Array[Pair[String, Pair[File, File]]] per_chr_bam_bai_and_id

        Boolean is_ont
        Boolean is_r10_4_pore_or_later
        String model_for_dv_andor_pepper

        # reference info
        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File? ref_scatter_interval_list_locator
        File? ref_scatter_interval_list_ids

        # smallVar-specific args
        Boolean run_clair3

        Int? dv_threads
        Int? dv_memory
        Boolean use_gpu = false

        # optimization
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    output {
        File? clair_vcf = MergeAndSortClairVCFs.vcf
        File? clair_tbi = MergeAndSortClairVCFs.tbi
        File? clair_gvcf = MergeAndSortClair_gVCFs.vcf
        File? clair_gtbi = MergeAndSortClair_gVCFs.tbi

        File haplotagged_bam = use_this_haptag_bam
        File haplotagged_bai = use_this_haptag_bai

        File dv_g_vcf = use_this_dv_g_vcf
        File dv_g_tbi = use_this_dv_g_tbi

        File dv_phased_vcf = use_this_dv_phased_vcf
        File dv_phased_tbi = use_this_dv_phased_tbi

        File dv_vcf_phasing_stats_tsv = use_this_phasing_stats_tsv
        File dv_vcf_phasing_stats_gtf = use_this_phasing_stats_gtf
    }

    #################################
    # shard WG specifically for small variant calling
    #################################
    # but if custom sharding isn't requested, then per-chr sharding is already done, so no need to to redo
    if (defined(ref_scatter_interval_list_locator)) {
        call ShardWholeGenome.Split as CustomSplitBamForSmallVar {
            input:
                ref_dict = ref_dict,
                bam = bam,
                bai = bai,
                ref_scatter_interval_list_locator = ref_scatter_interval_list_locator,
                ref_scatter_interval_list_ids = ref_scatter_interval_list_ids
        }
    }
    Array[Pair[String, Pair[File, File]]] how_to_shard_wg_for_calling = select_first([CustomSplitBamForSmallVar.id_bam_bai_of_shards,
                                                                                      per_chr_bam_bai_and_id])

    #################################
    # run clair3, if so requested
    #################################
    if (run_clair3) {
        scatter (triplet in how_to_shard_wg_for_calling) {
            call Clair3.Clair {
                input:
                    bam = triplet.right.left,
                    bai = triplet.right.right,

                    ref_fasta     = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,

                    preset = if is_ont then "ONT" else "CCS",
                    zones = zones
            }
        }
        call VariantUtils.MergeAndSortVCFs as MergeAndSortClairVCFs {
            input:
                vcfs = Clair.vcf,
                ref_fasta_fai = ref_fasta_fai,
                prefix = prefix + ".clair"
        }
        call VariantUtils.MergeAndSortVCFs as MergeAndSortClair_gVCFs {
            input:
                vcfs = Clair.gvcf,
                ref_fasta_fai = ref_fasta_fai,
                prefix = prefix + ".clair.g"
        }
    }

    #################################
    # DV, major workhorse (except for legacy ONT data)
    #################################
    if ((!is_ont) || is_r10_4_pore_or_later) { # pacbio or recent ONT data

        #############
        # per shard DV
        scatter (triplet in how_to_shard_wg_for_calling) {
            File shard_bam = triplet.right.left
            File shard_bai = triplet.right.right

            if (!use_gpu) {
                call DeepVariant.DV as DeepV {
                    input:
                        bam           = shard_bam,
                        bai           = shard_bai,
                        ref_fasta     = ref_fasta,
                        ref_fasta_fai = ref_fasta_fai,

                        model_type = model_for_dv_andor_pepper,

                        threads = select_first([dv_threads]),
                        memory  = select_first([dv_memory]),
                        zones = zones
                }
            }

            if (use_gpu) {
                call DeepVariant.DV_gpu as DeepV_G {
                    input:
                        bam           = shard_bam,
                        bai           = shard_bai,
                        ref_fasta     = ref_fasta,
                        ref_fasta_fai = ref_fasta_fai,

                        model_type = model_for_dv_andor_pepper,

                        threads = select_first([dv_threads]),
                        memory  = select_first([dv_memory]),
                        zones = zones
                }
            }

            File dv_vcf = select_first([DeepV.VCF, DeepV_G.VCF])
            File dv_gvcf = select_first([DeepV.gVCF, DeepV_G.gVCF])
        }

        String dv_prefix = prefix + ".deepvariant"

        #############
        # merge
        call VariantUtils.MergeAndSortVCFs as MergeDeepVariantGVCFs {
            input:
                vcfs     = dv_gvcf,
                prefix   = dv_prefix + ".g",
                ref_fasta_fai = ref_fasta_fai
        }

        call VariantUtils.MergeAndSortVCFs as MergeDeepVariantVCFs {
            input:
                vcfs     = dv_vcf,
                prefix   = dv_prefix,
                ref_fasta_fai = ref_fasta_fai
        }

        call Margin.MarginPhase as MarginPhase {
            input:
                bam=bam, bai=bai, ref_fasta=ref_fasta, ref_fasta_fai=ref_fasta_fai,
                unphased_vcf=MergeDeepVariantVCFs.vcf, unphased_tbi=MergeDeepVariantVCFs.tbi,
                zones = zones,
                data_type= if (is_ont) then "ONT" else "PacBio"
        }
        call WhatsHap.Stats as MPhaseStats  { input: phased_vcf=MarginPhase.phased_vcf, phased_tbi=MarginPhase.phased_tbi}

        #############
        # phase variant using whatshap phase, because it is much slower than margin phase, we do scatter-gather
        scatter (triplet in per_chr_bam_bai_and_id) {
            call VariantUtils.SubsetVCF as ChopDVVCF {
                input:
                    vcf_gz  = MergeDeepVariantVCFs.vcf,
                    vcf_tbi = MergeDeepVariantVCFs.tbi,
                    locus = triplet.left,
                    prefix = basename(MergeDeepVariantVCFs.vcf, ".vcf.gz") + "." + triplet.left
            }
            call WhatsHap.Phase as WhatsHapPhase {
                input :
                    chromosome=triplet.left,
                    bam=triplet.right.left, bai=triplet.right.right,
                    unphased_vcf=ChopDVVCF.subset_vcf, unphased_tbi=ChopDVVCF.subset_tbi,
                    ref_fasta=ref_fasta, ref_fasta_fai=ref_fasta_fai,
            }
        }
        call VariantUtils.MergePerChrCalls as MergeWhatsHapPhasedPerChrVCFs {
            input:
            vcfs = WhatsHapPhase.phased_vcf, ref_dict = ref_dict, prefix = basename(MergeDeepVariantVCFs.vcf, ".vcf.gz") + ".whatshap-phased"
        }
        call WhatsHap.Stats as WHPhaseStats { input: phased_vcf=MergeWhatsHapPhasedPerChrVCFs.vcf, phased_tbi=MergeWhatsHapPhasedPerChrVCFs.tbi}

        #############
        # if we later decide whatshap is better, change these
        File phased_snp_vcf = MarginPhase.phased_vcf
        File phased_snp_tbi = MarginPhase.phased_tbi
        File phasing_stats_tsv = MPhaseStats.stats_tsv
        File phasing_stats_gtf = MPhaseStats.stats_gtf

        #############
        # haplotag
        call WhatsHap.HaploTagBam as WhatsHapTag {
            input:
                to_tag_bam = bam,
                to_tag_bai = bai,
                ref_fasta  = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                phased_vcf = phased_snp_vcf,
                phased_tbi = phased_snp_tbi
        }
    }

    #################################
    # PEPPER-Margin-DeepVariant is the ONT data is before the R10.4 simplex/duplex generation
    #################################
    # old ONT data
    if (is_ont && (! is_r10_4_pore_or_later)) {
        scatter (triplet in how_to_shard_wg_for_calling) {
            call ONTPepper.Pepper {
                input:
                    bam           = triplet.right.left,
                    bai           = triplet.right.right,
                    ref_fasta     = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,
                    model = model_for_dv_andor_pepper,
                    threads       = select_first([dv_threads]),
                    memory        = select_first([dv_memory]),
                    zones = zones
            }
        }

        String pepper_prefix = prefix + ".PEPPER_deepvariant_margin"

        call VariantUtils.MergeAndSortVCFs as MergePEPPERGVCFs {
            input:
                vcfs     = Pepper.gVCF,
                prefix   = pepper_prefix + ".g",
                ref_fasta_fai = ref_fasta_fai
        }

        call VariantUtils.MergeAndSortVCFs as MergePEPPERPhasedVCFs {
            input:
                vcfs     = Pepper.phasedVCF,
                prefix   = pepper_prefix + ".phased",
                ref_fasta_fai = ref_fasta_fai
        }

        call Utils.MergeBams as MergePEPPERHapTaggedBam {
            input:
                bams = Pepper.hap_tagged_bam,
                prefix = prefix +  ".MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged"
        }

        call WhatsHap.Stats as ONTPhaseStatsLegacy  { input: phased_vcf=MergePEPPERPhasedVCFs.vcf, phased_tbi=MergePEPPERPhasedVCFs.tbi}

        File legacy_ont_dvp_g_vcf = MergePEPPERGVCFs.vcf
        File legacy_ont_dvp_g_tbi = MergePEPPERGVCFs.tbi
        File legacy_ont_dvp_phased_vcf = MergePEPPERPhasedVCFs.vcf
        File legacy_ont_dvp_phased_tbi = MergePEPPERPhasedVCFs.tbi
        File legacy_ont_dvp_haplotagged_bam = MergePEPPERHapTaggedBam.merged_bam
        File legacy_ont_dvp_haplotagged_bai = MergePEPPERHapTaggedBam.merged_bai
        File legacy_ont_dvp_phased_vcf_stats_tsv = ONTPhaseStatsLegacy.stats_tsv
        File legacy_ont_dvp_phased_vcf_stats_gtf = ONTPhaseStatsLegacy.stats_gtf
    }

    File use_this_dv_g_vcf = select_first([MergeDeepVariantGVCFs.vcf, legacy_ont_dvp_g_vcf])
    File use_this_dv_g_tbi = select_first([MergeDeepVariantGVCFs.tbi, legacy_ont_dvp_g_tbi])

    File use_this_dv_phased_vcf = select_first([phased_snp_vcf, legacy_ont_dvp_phased_vcf])
    File use_this_dv_phased_tbi = select_first([phased_snp_tbi, legacy_ont_dvp_phased_tbi])

    File use_this_phasing_stats_tsv = select_first([phasing_stats_tsv, legacy_ont_dvp_phased_vcf_stats_tsv])
    File use_this_phasing_stats_gtf = select_first([phasing_stats_gtf, legacy_ont_dvp_phased_vcf_stats_gtf])

    File use_this_haptag_bam = select_first([WhatsHapTag.tagged_bam, legacy_ont_dvp_haplotagged_bam])
    File use_this_haptag_bai = select_first([WhatsHapTag.tagged_bai, legacy_ont_dvp_haplotagged_bai])
}
