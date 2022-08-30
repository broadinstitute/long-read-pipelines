version 1.0

import "Utils.wdl"
import "VariantUtils.wdl"

import "PBSV.wdl"
import "Sniffles2.wdl" as Sniffles2

import "Clair.wdl" as Clair3

import "CCSPepper.wdl"

workflow CallVariants {
    meta {
        descrition: "A workflow for calling small and/or structural variants from an aligned CCS BAM file."
    }
    input {
        File bam
        File bai
        Int minsvlen
        String prefix
        String sample_id


        File ref_fasta
        File ref_fasta_fai
        File ref_dict


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
            String contig_for_small_var = c[0]

            call Utils.SubsetBam as SmallVariantsScatter {
                input:
                    bam = bam,
                    bai = bai,
                    locus = contig_for_small_var
            }

            call Clair3.Clair {
                input:
                    bam = SmallVariantsScatter.subset_bam,
                    bai = SmallVariantsScatter.subset_bai,

                    ref_fasta     = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,

                    sites_vcf = sites_vcf,
                    sites_vcf_tbi = sites_vcf_tbi,

                    preset = "CCS",
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
        # todo: phasing isn't done for CCS data yet, waiting for Pepper Team to respond
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

                call CCSPepper.CCSPepper {
                    input:
                        bam           = size_balanced_scatter.subset_bam,
                        bai           = size_balanced_scatter.subset_bai,
                        ref_fasta     = ref_fasta,
                        ref_fasta_fai = ref_fasta_fai,

                        pepper_threads = select_first([dvp_threads]),
                        pepper_memory  = select_first([dvp_memory]),
                        dv_threads = select_first([dvp_threads]),
                        dv_memory  = select_first([dvp_memory]),
                        zones = arbitrary.zones
                }
            }

            String dvp_prefix = prefix + ".deepvariant_pepper"

            call VariantUtils.MergeAndSortVCFs as MergeDeepVariantGVCFs {
                input:
                    vcfs     = CCSPepper.gVCF,
                    prefix   = dvp_prefix + ".g",
                    ref_fasta_fai = ref_fasta_fai
            }

            # todo: phasing VCF could happen here, i.e. on gathered VCFs as that's going to be less intensive
            call VariantUtils.MergeAndSortVCFs as MergeDeepVariantVCFs {
                input:
                    vcfs     = CCSPepper.VCF,
                    prefix   = dvp_prefix,
                    ref_fasta_fai = ref_fasta_fai
            }

            call Utils.MergeBams {
                input:
                    bams = CCSPepper.hap_tagged_bam,
                    prefix = prefix +  ".MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged"
            }

            call CCSPepper.MarginPhase {
                input:
                    bam           = bam,
                    bai           = bai,
                    unphased_vcf  = MergeDeepVariantVCFs.vcf,
                    unphased_vcf_tbi = MergeDeepVariantVCFs.tbi,
                    ref_fasta     = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,
                    memory        = select_first([dvp_memory, 64]),
                    zones = arbitrary.zones
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
                String contig_for_sv = c[0]

                call Utils.SubsetBam {
                    input:
                        bam = bam,
                        bai = bai,
                        locus = contig_for_sv
                }

                call PBSV.RunPBSV {
                    input:
                        bam = SubsetBam.subset_bam,
                        bai = SubsetBam.subset_bai,
                        ref_fasta = ref_fasta,
                        ref_fasta_fai = ref_fasta_fai,
                        prefix = prefix,
                        tandem_repeat_bed = tandem_repeat_bed,
                        is_ccs = true,
                        zones = arbitrary.zones
                }

                call Sniffles2.sample_sv as sample_snf_fast {
                    input:
                        bam    = SubsetBam.subset_bam,
                        bai    = SubsetBam.subset_bai,
                        minsvlen = minsvlen,
                        chr    = contig_for_sv,
                        sample_id = sample_id,
                        prefix = prefix
                }

                call Utils.InferSampleName {
                    input:
                        bam = SubsetBam.subset_bam,
                        bai = SubsetBam.subset_bai
                }
#                call VariantUtils.FixSnifflesVCF {
#                    input:
#                        vcf = Sniffles.vcf,
#                        sample_name = InferSampleName.sample_name
#                }
            }

            call VariantUtils.MergePerChrCalls as MergePBSVVCFs {
                input:
                    vcfs     = RunPBSV.vcf,
                    ref_dict = ref_dict,
                    prefix   = prefix + ".pbsv"
            }

#            call VariantUtils.CollectDefinitions as UnionHeadersSnifflesVCFs {
#                input:
#                    vcfs = FixSnifflesVCF.sortedVCF
#            }
#            call VariantUtils.MergeAndSortVCFs as MergeSnifflesVCFs {
#                input:
#                    vcfs   = FixSnifflesVCF.sortedVCF,
#                    ref_fasta_fai = ref_fasta_fai,
#                    header_definitions_file = UnionHeadersSnifflesVCFs.union_definitions,
#                    prefix = prefix + ".sniffles"
#            }
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
                    is_ccs = true,
                    zones = arbitrary.zones
            }
            call VariantUtils.ZipAndIndexVCF as ZipAndIndexPBSV {input: vcf = PBSVslow.vcf }

            call Sniffles2.sample_sv as sample_snf_slow {
                input:
                    bam    = bam,
                    bai    = bai,
                    minsvlen = minsvlen,
                    sample_id = sample_id,
                    prefix = prefix
            }


#            call Utils.InferSampleName as infer {input: bam = bam, bai = bai}
#            call VariantUtils.FixSnifflesVCF as ZipAndIndexSniffles {input: vcf = SnifflesSlow.vcf, sample_name = infer.sample_name}
        }
    }

    output {
        File? sniffles_vcf = select_first([sample_snf_fast.vcf, sample_snf_slow.vcf])
#        File? sniffles_tbi = select_first([MergeSnifflesVCFs.tbi, ZipAndIndexSniffles.tbi])
        File? sniffles_snf = select_first([sample_snf_fast.snf,sample_snf_slow.snf])

        File? pbsv_vcf = select_first([MergePBSVVCFs.vcf, ZipAndIndexPBSV.vcfgz])
        File? pbsv_tbi = select_first([MergePBSVVCFs.tbi, ZipAndIndexPBSV.tbi])

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
    }
}
