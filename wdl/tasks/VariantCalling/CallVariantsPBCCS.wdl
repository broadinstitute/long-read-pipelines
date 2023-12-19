version 1.0

import "../Utility/Utils.wdl"
import "../Utility/VariantUtils.wdl"
import "PBSV.wdl"
import "Sniffles2.wdl" as Sniffles2
import "Clair.wdl" as Clair3
import "CCSPepper.wdl" as Pepper

workflow CallVariants {

    meta {
        description: "A workflow for calling small and/or structural variants from an aligned CCS BAM file."
    }

    parameter_meta {
        bam: "Aligned CCS BAM file"
        bai: "Index for the aligned CCS BAM file"
        minsvlen: "Minimum SV length in bp (default: 50)"
        prefix: "Prefix for output file names"
        output_bucket: "Cloud path for output storage"
        sample_id: "Sample ID"
        ref_fasta: "Reference FASTA file"
        ref_fasta_fai: "Index for the reference FASTA file"
        ref_dict: "Dictionary for the reference FASTA file"
        regions_file: "TSV file of genomic regions to process (regions on each line processed together)"
        call_svs: "Call structural variants or not"
        fast_less_sensitive_sv: "to trade less sensitive SV calling for faster speed"
        tandem_repeat_bed: "BED file containing TRF finder for better PBSV calls (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"
        call_small_variants: "Call small variants or not"
        sites_vcf: "for use with Clair"
        sites_vcf_tbi: "for use with Clair"
        run_dv_pepper_analysis: "to turn on DV-Pepper analysis or not (non-trivial increase in cost and runtime)"
    }

    input {
        File bam
        File bai
        Int minsvlen = 50
        String prefix
        String? output_bucket
        String sample_id

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File regions_file

        Boolean call_svs
        Boolean fast_less_sensitive_sv
        File? tandem_repeat_bed

        Boolean call_small_variants
        File? sites_vcf
        File? sites_vcf_tbi

        Boolean run_dv_pepper_analysis
    }

    ######################################################################
    # Block for small variants handling
    ######################################################################

    # read a tsv file listing the chromosomal regions to process
    Array[Array[String]] regionsArray = read_tsv(regions_file)

    if (call_small_variants) {
        String? snp_dir = if (defined(output_bucket)) then sub(select_first([output_bucket]), "/$", "") + "/variants/small" else output_bucket
        scatter (regions in regionsArray) {
            call Clair3.Clair {
                input:
                    bam = bam,
                    bai = bai,
                    prefix = prefix,
                    ref_fasta = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,
                    sites_vcf = sites_vcf,
                    sites_vcf_tbi = sites_vcf_tbi,
                    regions = regions,
                    preset = "CCS"
            }
        }

        call VariantUtils.MergeAndSortVCFs as MergeAndSortClairVCFs {
            input:
                vcfs = Clair.vcf,
                ref_fasta_fai = ref_fasta_fai,
                prefix = prefix + ".clair",
                output_bucket = snp_dir
        }

        call VariantUtils.MergeAndSortVCFs as MergeAndSortClair_gVCFs {
            input:
                vcfs = Clair.gvcf,
                ref_fasta_fai = ref_fasta_fai,
                prefix = prefix + ".clair.g",
                output_bucket = snp_dir
        }

        # todo: phasing isn't done for CCS data yet, waiting for Pepper Team to respond
        if (run_dv_pepper_analysis) {
            scatter (regions in regionsArray) {
                call Pepper.CCSPepper {
                    input:
                        bam = bam,
                        bai = bai,
                        ref_fasta = ref_fasta,
                        ref_fasta_fai = ref_fasta_fai,
                        regions = regions
                }
            }

            String dvp_prefix = prefix + ".deepvariant_pepper"

            call VariantUtils.MergeAndSortVCFs as MergeDeepVariantGVCFs {
                input:
                    vcfs = CCSPepper.gVCF,
                    ref_fasta_fai = ref_fasta_fai,
                    prefix = dvp_prefix + ".g",
                    output_bucket = snp_dir
            }

            # todo: phasing VCF could happen here, i.e. on gathered VCFs as that's going to be less intensive
            call VariantUtils.MergeAndSortVCFs as MergeDeepVariantVCFs {
                input:
                    vcfs = CCSPepper.VCF,
                    ref_fasta_fai = ref_fasta_fai,
                    prefix = dvp_prefix,
                    output_bucket = snp_dir
            }

            call Utils.MergeBams {
                input:
                    bams = CCSPepper.hap_tagged_bam,
                    outputBamName = "~{prefix}.MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam",
                    outputBucket = snp_dir
            }

            call Pepper.MarginPhase {
                input:
                    bam = bam,
                    bai = bai,
                    unphased_vcf = MergeDeepVariantVCFs.vcf,
                    unphased_vcf_tbi = MergeDeepVariantVCFs.tbi,
                    ref_fasta = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,
                    output_bucket = snp_dir
            }
        }
    }

    ######################################################################
    # Block for SV handling
    ######################################################################
    if (call_svs) {
        String? sv_dir = if (defined(output_bucket)) then sub(select_first([output_bucket]), "/$", "" + "/variants/sv" else output_bucket
        if (fast_less_sensitive_sv) {
            scatter (regions in regionsArray) {
                call PBSV.RunPBSV {
                    input:
                        bam = bam,
                        bai = bai,
                        is_ccs = true,
                        ref_fasta = ref_fasta,
                        ref_fasta_fai = ref_fasta_fai,
                        prefix = prefix,
                        tandem_repeat_bed = tandem_repeat_bed,
                        regions = regions,
                        n_tasks = length(regionsArray)
                }

            }

            call VariantUtils.MergeAndSortVCFs as MergePBSVVCFs {
                input:
                    vcfs = RunPBSV.vcf,
                    ref_fasta_fai = ref_fasta_fai,
                    prefix = prefix + ".pbsv",
                    output_bucket = sv_dir
            }

        }

        if (!fast_less_sensitive_sv) {
            call PBSV.RunPBSV as PBSVslow {
                input:
                    bam = bam,
                    bai = bai,
                    is_ccs = true,
                    ref_fasta = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,
                    prefix = prefix,
                    output_bucket = sv_dir,
                    tandem_repeat_bed = tandem_repeat_bed
            }
        }

        call Sniffles2.SampleSV as Sniffles2SV {
            input:
                bam    = bam,
                bai    = bai,
                minsvlen = minsvlen,
                sample_id = sample_id,
                prefix = prefix
        }
    }

    output {
        File? sniffles_vcf = Sniffles2SV.vcf
        File? sniffles_tbi = Sniffles2SV.tbi
        File? sniffles_snf = Sniffles2SV.snf

        # can't do a select_first here because if call_svs is false, then both will be null and select_first will fail
        File? pbsv_vcf = if defined(MergePBSVVCFs.vcf) then MergePBSVVCFs.vcf else PBSVslow.vcf
        File? pbsv_tbi = if defined(MergePBSVVCFs.tbi) then MergePBSVVCFs.tbi else PBSVslow.tbi

        File? clair_vcf = MergeAndSortClairVCFs.vcf
        File? clair_tbi = MergeAndSortClairVCFs.tbi

        File? clair_gvcf = MergeAndSortClair_gVCFs.vcf
        File? clair_gtbi = MergeAndSortClair_gVCFs.tbi

        File? dvp_g_vcf = MergeDeepVariantGVCFs.vcf
        File? dvp_g_tbi = MergeDeepVariantGVCFs.tbi
        File? dvp_vcf = MergeDeepVariantVCFs.vcf
        File? dvp_tbi = MergeDeepVariantVCFs.tbi
        File? dvp_phased_vcf = MarginPhase.phasedVCF
        File? dvp_phased_tbi = MarginPhase.phasedtbi

        File? haplotagged_bam = MergeBams.merged_bam
    }
}
