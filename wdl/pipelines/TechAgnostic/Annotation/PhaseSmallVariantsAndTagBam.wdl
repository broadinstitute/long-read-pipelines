version 1.0

import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../../tasks/Utility/VariantUtils.wdl"

import "../../../tasks/Alignment/WhatsHap.wdl"
import "../../../tasks/VariantCalling/MarginPhase.wdl" as Margin

workflow Run {
    meta {
        desciption:
        "For read-based small variant VCF phasing and haplotagging BAM"
    }
    parameter_meta {
        use_margin_for_tagging: "if false, will use margin-phased VCF for haplotagging the BAM"
    }

    input {
        Boolean use_margin_for_tagging

        File bam
        File bai
        Array[Pair[String, Pair[File, File]]] per_chr_bam_bai_and_id
        Boolean is_ont

        File unphased_vcf
        File unphased_tbi

        Map[String, String] ref_map

        String gcs_variants_out_dir
        String gcs_tagged_bam_out_dir

        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }

    output {
        File margin_phased_vcf = FinalizeDVMarginPhasedVcf.gcs_path
        File margin_phased_tbi = FinalizeDVMarginPhasedTbi.gcs_path
        File margin_phasing_stats_tsv = FinalizeDVMarginPhasedVcfStatusTSV.gcs_path
        File margin_phasing_stats_gtf = FinalizeDVMarginPhasedVcfStatusGtf.gcs_path

        File whatshap_phased_vcf = FinalizeDVWhatsHapPhasedVcf.gcs_path
        File whatshap_phased_tbi = FinalizeDVWhatsHapPhasedTbi.gcs_path
        File whatshap_phasing_stats_tsv = FinalizeDVWhatsHapPhasedVcfStatusTSV.gcs_path
        File whatshap_phasing_stats_gtf = FinalizeDVWhatsHapPhasedVcfStatusGtf.gcs_path

        File hap_tagged_bam = FinalizeHapTaggedBam.gcs_path
        File hap_tagged_bai = FinalizeHapTaggedBai.gcs_path
        String haplotagged_bam_tagger = if use_margin_for_tagging then "MARGIN" else "WhatsHap"
    }

    ####################################################################################################################################
    # TODO: we need to comeback and make a choice on which phasing tool to use
    ####################################################################################################################################
    #############
    # phase variants using margin
    #############
    call Margin.MarginPhase as MarginPhase {
        input:
            bam=bam, bai=bai,
            data_type= if (is_ont) then "ONT" else "PacBio",

            unphased_vcf=unphased_vcf,
            unphased_tbi=unphased_tbi,

            ref_fasta=ref_map['fasta'], ref_fasta_fai=ref_map['fai'],
            zones = zones
    }
    call WhatsHap.Stats as MPhaseStats  { input: phased_vcf=MarginPhase.phased_vcf, phased_tbi=MarginPhase.phased_tbi}

    call FF.FinalizeToFile as FinalizeDVMarginPhasedVcf          { input: outdir = gcs_variants_out_dir, file = MarginPhase.phased_vcf }
    call FF.FinalizeToFile as FinalizeDVMarginPhasedTbi          { input: outdir = gcs_variants_out_dir, file = MarginPhase.phased_tbi }
    call FF.FinalizeToFile as FinalizeDVMarginPhasedVcfStatusTSV { input: outdir = gcs_variants_out_dir, file = MPhaseStats.stats_tsv }
    call FF.FinalizeToFile as FinalizeDVMarginPhasedVcfStatusGtf { input: outdir = gcs_variants_out_dir, file = MPhaseStats.stats_gtf }

    #############
    # phase variant using whatshap phase; but because it is much slower than margin phase, we do scatter-gather
    #############
    scatter (triplet in per_chr_bam_bai_and_id) {
        call VariantUtils.SubsetVCF as ChopDVVCF {
            input:
                vcf_gz  = unphased_vcf,
                vcf_tbi = unphased_tbi,
                locus = triplet.left,
                prefix = basename(unphased_vcf, ".vcf.gz") + "." + triplet.left
        }
        call WhatsHap.Phase as WhatsHapPhase {
            input :
                chromosome=triplet.left,
                bam=triplet.right.left, bai=triplet.right.right,
                unphased_vcf=ChopDVVCF.subset_vcf, unphased_tbi=ChopDVVCF.subset_tbi,
                ref_fasta=ref_map['fasta'], ref_fasta_fai=ref_map['fai'],
        }
    }
    call VariantUtils.MergePerChrCalls as MergeWhatsHapPhasedPerChrVCFs {
        input:
        vcfs = WhatsHapPhase.phased_vcf, ref_dict = ref_map['dict'],
        prefix = basename(unphased_vcf, ".vcf.gz") + ".whatshap-phased"
    }
    call WhatsHap.Stats as WHPhaseStats { input: phased_vcf=MergeWhatsHapPhasedPerChrVCFs.vcf, phased_tbi=MergeWhatsHapPhasedPerChrVCFs.tbi}

    call FF.FinalizeToFile as FinalizeDVWhatsHapPhasedVcf          { input: outdir = gcs_variants_out_dir, file = MergeWhatsHapPhasedPerChrVCFs.vcf }
    call FF.FinalizeToFile as FinalizeDVWhatsHapPhasedTbi          { input: outdir = gcs_variants_out_dir, file = MergeWhatsHapPhasedPerChrVCFs.tbi }
    call FF.FinalizeToFile as FinalizeDVWhatsHapPhasedVcfStatusTSV { input: outdir = gcs_variants_out_dir, file = WHPhaseStats.stats_tsv }
    call FF.FinalizeToFile as FinalizeDVWhatsHapPhasedVcfStatusGtf { input: outdir = gcs_variants_out_dir, file = WHPhaseStats.stats_gtf }

    ####################################################################################################################################
    #############
    # !CHOICE! haplotag with WhatsHap, but using which phased VCF?!
    #############
    File phased_snp_vcf = if (use_margin_for_tagging) then MarginPhase.phased_vcf else MergeWhatsHapPhasedPerChrVCFs.vcf
    File phased_snp_tbi = if (use_margin_for_tagging) then MarginPhase.phased_tbi else MergeWhatsHapPhasedPerChrVCFs.tbi

    call WhatsHap.HaploTagBam as WhatsHapTag {
        input:
            to_tag_bam = bam,
            to_tag_bai = bai,
            ref_fasta  = ref_map['fasta'],
            ref_fasta_fai = ref_map['fai'],
            phased_vcf = phased_snp_vcf,
            phased_tbi = phased_snp_tbi
    }

    call FF.FinalizeToFile as FinalizeHapTaggedBam { input: outdir = gcs_tagged_bam_out_dir, file = WhatsHapTag.tagged_bam }
    call FF.FinalizeToFile as FinalizeHapTaggedBai { input: outdir = gcs_tagged_bam_out_dir, file = WhatsHapTag.tagged_bai }
}
