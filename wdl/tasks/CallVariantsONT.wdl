## #

version 1.0

import "Utils.wdl"
import "VariantUtils.wdl"

import "PBSV.wdl"
import "Sniffles.wdl"

import "Clair.wdl" as Clair3

import "ONTPepper.wdl"

workflow CallVariants {

    meta {
        description: "A workflow for calling small and/or structural variants from an aligned ONT BAM file."
    }

    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String prefix

        Boolean call_svs
        Boolean fast_less_sensitive_sv
        File? tandem_repeat_bed

        Boolean call_small_variants
        Boolean call_small_vars_on_mitochondria
        File? sites_vcf
        File? sites_vcf_tbi

        Boolean run_dv_pepper_analysis
        Int? dvp_threads
        Int? dvp_memory
        File? ref_scatter_interval_list_locator
        File? ref_scatter_interval_list_ids
    }

    parameter_meta {
        fast_less_sensitive_sv:  "to trade less sensitive SV calling for faster speed"
        tandem_repeat_bed:       "BED file containing TRF finder for better PBSV calls (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"

        call_small_vars_on_mitochondria: "if false, will not attempt to call variants on mitochondria"
        sites_vcf:     "for use with Clair"
        sites_vcf_tbi: "for use with Clair"

        run_dv_pepper_analysis:  "to turn on DV-Pepper analysis or not (non-trivial increase in cost and runtime)"
        ref_scatter_interval_list_locator: "A file holding paths to interval_list files; needed only when running DV-Pepper"
        ref_scatter_interval_list_ids:     "A file that gives short IDs to the interval_list files; needed only when running DV-Pepper"
    }

    ######################################################################
    # Block for small variants handling
    ######################################################################

    call Utils.RandomZoneSpewer as arbitrary {input: num_of_zones = 3}

    # todo: merge the two scattering scheme into a better one
    if (call_small_variants) {
        # Scatter by chromosome
        Array[String] default_filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']
        Array[String] use_filter = if (call_small_vars_on_mitochondria) then default_filter else flatten([['chrM'],default_filter])
        call Utils.MakeChrIntervalList as SmallVariantsScatterPrepp {
            input:
                ref_dict = ref_dict,
                filter = use_filter
        }

        scatter (c in SmallVariantsScatterPrepp.chrs) {
            String chr = c[0]

            call Utils.SubsetBam as SmallVariantsScatter {
                input:
                    bam = bam,
                    bai = bai,
                    locus = chr
            }

            call Clair3.Clair {
                input:
                    bam = SmallVariantsScatter.subset_bam,
                    bai = SmallVariantsScatter.subset_bai,

                    ref_fasta     = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,

                    sites_vcf = sites_vcf,
                    sites_vcf_tbi = sites_vcf_tbi,

                    preset = "ONT",
                    zones = arbitrary.zones
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

        # size-balanced scatter
        if (run_dv_pepper_analysis) {
            File scatter_interval_list_ids = select_first([ref_scatter_interval_list_ids])
            File scatter_interval_list_loc = select_first([ref_scatter_interval_list_locator])
            Array[String] interval_list_ids   = read_lines(scatter_interval_list_ids)
            Array[String] interval_list_files = read_lines(scatter_interval_list_loc)
            Array[Pair[String, String]] ided_interval_list_files = zip(interval_list_ids, interval_list_files)

            scatter (pair in ided_interval_list_files) {
                call Utils.ResilientSubsetBam as size_balanced_scatter {
                    input:
                        bam = bam,
                        bai = bai,
                        interval_list_file = pair.right,
                        interval_id = pair.left,
                        prefix = basename(bam, ".bam")
                }

                call ONTPepper.Pepper {
                    input:
                        bam           = size_balanced_scatter.subset_bam,
                        bai           = size_balanced_scatter.subset_bai,
                        ref_fasta     = ref_fasta,
                        ref_fasta_fai = ref_fasta_fai,
                        threads       = select_first([dvp_threads]),
                        memory        = select_first([dvp_memory]),
                        zones = arbitrary.zones
                }
            }

            String dvp_prefix = prefix + ".deepvariant_pepper"

            call VariantUtils.MergeAndSortVCFs as MergeDeepVariantGVCFs {
                input:
                    vcfs     = Pepper.gVCF,
                    prefix   = dvp_prefix + ".g",
                    ref_fasta_fai = ref_fasta_fai
            }

            call VariantUtils.MergeAndSortVCFs as MergeDeepVariantPhasedVCFs {
                input:
                    vcfs     = Pepper.phasedVCF,
                    prefix   = dvp_prefix + ".phased",
                    ref_fasta_fai = ref_fasta_fai
            }

            call VariantUtils.MergeAndSortVCFs as MergeDeepVariantVCFs {
                input:
                    vcfs     = Pepper.VCF,
                    prefix   = dvp_prefix,
                    ref_fasta_fai = ref_fasta_fai
            }

            call Utils.MergeBams {
                input:
                    bams = Pepper.hap_tagged_bam,
                    prefix = prefix +  ".MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged"
            }
        }
    }
    ######################################################################
    # Block for SV handling
    ######################################################################
    if (call_svs) {
        if (fast_less_sensitive_sv) {

            call Utils.MakeChrIntervalList {
            input:
                ref_dict = ref_dict,
                filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']
            }

            scatter (c in MakeChrIntervalList.chrs) {
                String contig = c[0]

                call Utils.SubsetBam {
                    input:
                        bam = bam,
                        bai = bai,
                        locus = contig
                }

                call PBSV.RunPBSV {
                    input:
                        bam = SubsetBam.subset_bam,
                        bai = SubsetBam.subset_bai,
                        ref_fasta = ref_fasta,
                        ref_fasta_fai = ref_fasta_fai,
                        prefix = prefix,
                        tandem_repeat_bed = tandem_repeat_bed,
                        is_ccs = false,
                        zones = arbitrary.zones
                }

                call Sniffles.Sniffles {
                    input:
                        bam    = SubsetBam.subset_bam,
                        bai    = SubsetBam.subset_bai,
                        chr    = contig,
                        prefix = prefix
                }

                call Utils.InferSampleName {
                    input:
                        bam = SubsetBam.subset_bam,
                        bai = SubsetBam.subset_bai
                }
                call VariantUtils.FixSnifflesVCF {
                    input:
                        vcf = Sniffles.vcf,
                        sample_name = InferSampleName.sample_name
                }
            }

            call VariantUtils.MergePerChrCalls as MergePBSVVCFs {
                input:
                    vcfs     = RunPBSV.vcf,
                    ref_dict = ref_dict,
                    prefix   = prefix + ".pbsv"
            }

            call VariantUtils.CollectDefinitions as UnionHeadersSnifflesVCFs {
                input:
                    vcfs = FixSnifflesVCF.sortedVCF
            }
            call VariantUtils.MergeAndSortVCFs as MergeSnifflesVCFs {
                input:
                    vcfs   = FixSnifflesVCF.sortedVCF,
                    ref_fasta_fai = ref_fasta_fai,
                    header_definitions_file = UnionHeadersSnifflesVCFs.union_definitions,
                    prefix = prefix + ".sniffles"
            }
        }

        if (!fast_less_sensitive_sv) {

            call PBSV.RunPBSV as PBSVslow {
                input:
                    bam = bam,
                    bai = bai,
                    ref_fasta = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,
                    prefix = prefix,
                    tandem_repeat_bed = tandem_repeat_bed,
                    is_ccs = false,
                    zones = arbitrary.zones
            }
            call VariantUtils.ZipAndIndexVCF as ZipAndIndexPBSV {input: vcf = PBSVslow.vcf }

            call Sniffles.Sniffles as SnifflesSlow {
                input:
                    bam    = bam,
                    bai    = bai,
                    prefix = prefix
            }
            call Utils.InferSampleName as infer {input: bam = bam, bai = bai}
            call VariantUtils.FixSnifflesVCF as ZipAndIndexSniffles {input: vcf = SnifflesSlow.vcf, sample_name = infer.sample_name}
        }
    }

    output {
        File? sniffles_vcf = select_first([MergeSnifflesVCFs.vcf, ZipAndIndexSniffles.sortedVCF])
        File? sniffles_tbi = select_first([MergeSnifflesVCFs.tbi, ZipAndIndexSniffles.tbi])

        File? pbsv_vcf = select_first([MergePBSVVCFs.vcf, ZipAndIndexPBSV.vcfgz])
        File? pbsv_tbi = select_first([MergePBSVVCFs.tbi, ZipAndIndexPBSV.tbi])

        File? clair_vcf = MergeAndSortClairVCFs.vcf
        File? clair_tbi = MergeAndSortClairVCFs.tbi

        File? clair_gvcf = MergeAndSortClair_gVCFs.vcf
        File? clair_gtbi = MergeAndSortClair_gVCFs.tbi

        File? dvp_phased_vcf = MergeDeepVariantPhasedVCFs.vcf
        File? dvp_phased_tbi = MergeDeepVariantPhasedVCFs.tbi
        File? dvp_g_vcf = MergeDeepVariantGVCFs.vcf
        File? dvp_g_tbi = MergeDeepVariantGVCFs.tbi
        File? dvp_vcf = MergeDeepVariantVCFs.vcf
        File? dvp_tbi = MergeDeepVariantVCFs.tbi
    }
}
