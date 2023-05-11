version 1.0

import "../Utility/Utils.wdl"
import "../Utility/VariantUtils.wdl"

import "ShardWholeGenome.wdl"

import "PBSV.wdl"
import "Sniffles2.wdl" as Sniffles2

import "Clair.wdl" as Clair3

import "DeepVariant.wdl"
import "../Alignment/WhatsHap.wdl"
import "MarginPhase.wdl" as Margin

workflow CallVariants {

    meta {
        description: "A workflow for calling small and/or structural variants from an aligned CCS BAM file."
    }

    parameter_meta {
        bam: "Aligned CCS BAM file"
        bai: "Index for the aligned CCS BAM file"
        prefix: "Prefix for output files"
        sample_id: "Sample ID"
        ref_fasta: "Reference FASTA file"
        ref_fasta_fai: "Index for the reference FASTA file"
        ref_dict: "Dictionary for the reference FASTA file"

        call_svs: "Call structural variants or not"
        minsvlen: "Minimum SV length in bp (default: 50)"
        pbsv_discover_per_chr: "Run the discover stage of PBSV per chromosome"
        tandem_repeat_bed: "BED file containing TRF finder for better SV calls (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"

        call_small_variants: "Call small variants or not"

        run_clair3: "to turn on Clair3 analysis or not (non-trivial increase in cost and runtime)"

        dv_threads: "number of threads for DeepVariant"
        dv_memory:  "memory for DeepVariant"
        use_gpu: "to use GPU acceleration or not on DeepVariant"
        ref_scatter_interval_list_locator: "A file holding paths to interval_list files; when not provided, will shard WG by contig (possibly slower)"
        ref_scatter_interval_list_ids: "A file that gives short IDs to the interval_list files; when not provided, will shard WG by contig (possibly slower)"
    }

    input {
        # sample info
        File bam
        File bai
        String prefix
        String sample_id

        # reference info
        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File? ref_scatter_interval_list_locator
        File? ref_scatter_interval_list_ids

        # sv-specific args
        Boolean call_svs
        Int minsvlen = 50
        File? tandem_repeat_bed

        Boolean pbsv_discover_per_chr

        # smallVar-specific args
        Boolean call_small_variants

        Boolean run_clair3

        Int? dv_threads
        Int? dv_memory
        Boolean use_gpu = false

        # optimization
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }

    ######################################################################
    # Block for small variants handling
    ######################################################################
    if (call_small_variants) {

        call ShardWholeGenome.Split as SplitBamForSmallVar {
            input:
                ref_dict = ref_dict,
                bam = bam,
                bai = bai,
                ref_scatter_interval_list_locator = ref_scatter_interval_list_locator,
                ref_scatter_interval_list_ids = ref_scatter_interval_list_ids
        }

        #################################
        # run clair3, if so requested
        #################################
        if (run_clair3) {
            scatter (pair in SplitBamForSmallVar.sharded_bam_bais) {
                call Clair3.Clair {
                    input:
                        bam = pair.left,
                        bai = pair.right,

                        ref_fasta     = ref_fasta,
                        ref_fasta_fai = ref_fasta_fai,

                        preset = "CCS",
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
        # running DeepVariant
        #################################
        scatter (pair in SplitBamForSmallVar.sharded_bam_bais) {
            if (!use_gpu) {
                call DeepVariant.DV as DeepV {
                    input:
                        bam           = pair.left,
                        bai           = pair.right,
                        ref_fasta     = ref_fasta,
                        ref_fasta_fai = ref_fasta_fai,

                        model_type = "PACBIO",

                        threads = select_first([dv_threads]),
                        memory  = select_first([dv_memory]),
                        zones = zones
                }
            }

            if (use_gpu) {
                call DeepVariant.DV_gpu as DeepV_G {
                    input:
                        bam           = pair.left,
                        bai           = pair.right,
                        ref_fasta     = ref_fasta,
                        ref_fasta_fai = ref_fasta_fai,

                        model_type = "PACBIO",

                        threads = select_first([dv_threads]),
                        memory  = select_first([dv_memory]),
                        zones = zones
                }
            }
            File dv_vcf = select_first([DeepV.VCF, DeepV_G.VCF])
            File dv_gvcf = select_first([DeepV.gVCF, DeepV_G.gVCF])
        }

        String dv_prefix = prefix + ".deepvariant"

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

        #############
        # to-be-fixed
        # here we run two phasing methods, because right now DV doesn't seem to do what they claim yet: direct phasing
        call WhatsHap.Phase as WhatsHapPhase {
            input :
                bam=bam, bai=bai, ref_fasta=ref_fasta, ref_fasta_fai=ref_fasta_fai,
                unphased_vcf=MergeDeepVariantVCFs.vcf, unphased_tbi=MergeDeepVariantVCFs.tbi
        }
        call WhatsHap.Stats as WHPhaseStats { input: phased_vcf=WhatsHapPhase.phased_vcf, phased_tbi=WhatsHapPhase.phased_tbi}

        call Margin.MarginPhase as MarginPhase {
            input:
                bam=bam, bai=bai, ref_fasta=ref_fasta, ref_fasta_fai=ref_fasta_fai,
                unphased_vcf=MergeDeepVariantVCFs.vcf, unphased_tbi=MergeDeepVariantVCFs.tbi,
                zones = zones,
                data_type='PacBio'
        }
        call WhatsHap.Stats as MPhaseStats  { input: phased_vcf=MarginPhase.phased_vcf, phased_tbi=MarginPhase.phased_tbi}

        File phased_snp_vcf = MarginPhase.phased_vcf
        File phased_snp_tbi = MarginPhase.phased_tbi
        File phasing_stats_tsv = MPhaseStats.stats_tsv
        File phasing_stats_gtf = MPhaseStats.stats_gtf
        #############

        call WhatsHap.HaploTagBam {
            input:
                to_tag_bam = bam,
                to_tag_bai = bai,
                ref_fasta  = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                phased_vcf = phased_snp_vcf,
                phased_tbi = phased_snp_tbi
        }
    }

    ######################################################################
    # Block for SV handling
    ######################################################################
    if (call_svs) {
        if (pbsv_discover_per_chr) {

            call Utils.MakeChrIntervalList {
            input:
                ref_dict = ref_dict,
                filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']
            }

            scatter (c in MakeChrIntervalList.chrs) {
                String contig_for_sv = c[0]

                call Utils.SubsetBam {
                    input:
                        bam = bam,
                        bai = bai,
                        locus = contig_for_sv
                }

                call PBSV.Discover as pbsv_discover_chr {
                    input:
                        bam = SubsetBam.subset_bam,
                        bai = SubsetBam.subset_bai,
                        ref_fasta = ref_fasta,
                        ref_fasta_fai = ref_fasta_fai,
                        tandem_repeat_bed = tandem_repeat_bed,
                        chr = contig_for_sv,
                        prefix = prefix,
                        zones = zones
                }
            }

            call PBSV.Call as pbsv_wg_call {
                input:
                    svsigs = pbsv_discover_chr.svsig,
                    ref_fasta = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,
                    ccs = true,
                    prefix = prefix + ".pbsv",
                    zones = zones
            }

            call VariantUtils.ZipAndIndexVCF as ZipAndIndexFastPBSV {input: vcf = pbsv_wg_call.vcf }
        }

        if (!pbsv_discover_per_chr) {
            call PBSV.RunPBSV as PBSVslow {
                input:
                    bam = bam,
                    bai = bai,
                    ref_fasta = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,
                    prefix = prefix,
                    tandem_repeat_bed = tandem_repeat_bed,
                    is_ccs = true,
                    zones = zones
            }

            call VariantUtils.ZipAndIndexVCF as ZipAndIndexPBSV {input: vcf = PBSVslow.vcf }
        }

        call Sniffles2.SampleSV as Sniffles2SV {
            input:
                bam    = bam,
                bai    = bai,
                minsvlen = minsvlen,
                sample_id = sample_id,
                prefix = prefix,
                tandem_repeat_bed = tandem_repeat_bed
        }
    }

    output {
        File? sniffles_vcf = Sniffles2SV.vcf
        File? sniffles_tbi = Sniffles2SV.tbi
        File? sniffles_snf = Sniffles2SV.snf

        File? pbsv_vcf = select_first([ZipAndIndexFastPBSV.vcfgz, ZipAndIndexPBSV.vcfgz])
        File? pbsv_tbi = select_first([ZipAndIndexFastPBSV.tbi, ZipAndIndexPBSV.tbi])

        File? clair_vcf = MergeAndSortClairVCFs.vcf
        File? clair_tbi = MergeAndSortClairVCFs.tbi
        File? clair_gvcf = MergeAndSortClair_gVCFs.vcf
        File? clair_gtbi = MergeAndSortClair_gVCFs.tbi

        File? dv_g_vcf = MergeDeepVariantGVCFs.vcf
        File? dv_g_tbi = MergeDeepVariantGVCFs.tbi

        File? dv_phased_vcf = phased_snp_vcf
        File? dv_phased_tbi = phased_snp_tbi

        File? dv_vcf_phasing_stats_tsv = phasing_stats_tsv
        File? dv_vcf_phasing_stats_gtf = phasing_stats_gtf

        File? haplotagged_bam = HaploTagBam.tagged_bam
        File? haplotagged_bai = HaploTagBam.tagged_bai
    }
}
