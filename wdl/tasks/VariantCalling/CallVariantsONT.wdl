version 1.0

import "../Utility/Utils.wdl"
import "../Utility/VariantUtils.wdl"
import "PBSV.wdl"
import "Sniffles2.wdl" as Sniffles2
import "Clair.wdl" as Clair3
import "ONTPepper.wdl"

workflow CallVariants {

    meta {
        descrition: "A workflow for calling small and/or structural variants from an aligned ONT BAM file."
    }

    parameter_meta {
        bam: "ONT BAM file"
        bai: "ONT BAM index file"
        minsvlen: "Minimum SV length in bp (default: 50)"
        prefix: "Prefix for output files"
        sample_id: "Sample ID"

        ref_map_file:       "Table indicating reference sequence, auxillary file locations, and metadata."

        call_svs: "Call structural variants"
        fast_less_sensitive_sv: "to trade less sensitive SV calling for faster speed"
        tandem_repeat_bed: "BED file containing TRF finder for better PBSV calls (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"

        call_small_variants: "Call small variants"
        call_small_vars_on_mitochondria: "if false, will not attempt to call variants on mitochondria"
        sites_vcf: "for use with Clair"
        sites_vcf_tbi: "for use with Clair"

        run_dv_pepper_analysis: "to turn on DV-Pepper analysis or not (non-trivial increase in cost and runtime)"
        dvp_threads: "number of threads for DV-Pepper"
        dvp_memory: "memory for DV-Pepper"
        ref_scatter_interval_list_locator: "A file holding paths to interval_list files; needed only when running DV-Pepper"
        ref_scatter_interval_list_ids: "A file that gives short IDs to the interval_list files; needed only when running DV-Pepper"
    }

    input {
        File bam
        File bai
        Int minsvlen = 50
        String prefix
        String sample_id

        File ref_map_file

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

    # Read ref map into map data type so we can access its fields:
    Map[String, String] ref_map = read_map(ref_map_file)

    ######################################################################
    # Block for small variants handling
    ######################################################################

    call Utils.RandomZoneSpewer as PickGcpZones {input: num_of_zones = 3}

    # todo: merge the two scattering scheme into a better one
    if (call_small_variants) {
        # Scatter by chromosome
        Array[String] default_filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']
        Array[String] use_filter = if (call_small_vars_on_mitochondria) then default_filter else flatten([['chrM', ref_map["mt_chr_name"]], default_filter])
        call Utils.MakeChrIntervalList as SmallVariantsScatterPrepp {
            input:
                ref_dict = ref_map["dict"],
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

                    ref_fasta     = ref_map["fasta"],
                    ref_fasta_fai = ref_map["fai"],

                    sites_vcf = sites_vcf,
                    sites_vcf_tbi = sites_vcf_tbi,

                    preset = "ONT",
                    zones = PickGcpZones.zones
            }
        }

        call VariantUtils.MergeAndSortVCFs as MergeAndSortClairVCFs {
            input:
                vcfs = Clair.vcf,
                ref_fasta_fai = ref_map["fai"],
                prefix = prefix + ".clair"
        }

        call VariantUtils.MergeAndSortVCFs as MergeAndSortClair_gVCFs {
            input:
                vcfs = Clair.gvcf,
                ref_fasta_fai = ref_map["fai"],
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
                        ref_fasta     = ref_map["fasta"],
                        ref_fasta_fai = ref_map["fai"],
                        threads       = select_first([dvp_threads]),
                        memory        = select_first([dvp_memory]),
                        zones = PickGcpZones.zones
                }
            }

            String dvp_prefix = prefix + ".deepvariant_pepper"

            call VariantUtils.MergeAndSortVCFs as MergeDeepVariantGVCFs {
                input:
                    vcfs     = Pepper.gVCF,
                    prefix   = dvp_prefix + ".g",
                    ref_fasta_fai = ref_map["fai"]
            }

            call VariantUtils.MergeAndSortVCFs as MergeDeepVariantPhasedVCFs {
                input:
                    vcfs     = Pepper.phasedVCF,
                    prefix   = dvp_prefix + ".phased",
                    ref_fasta_fai = ref_map["fai"]
            }

            call VariantUtils.MergeAndSortVCFs as MergeDeepVariantVCFs {
                input:
                    vcfs     = Pepper.VCF,
                    prefix   = dvp_prefix,
                    ref_fasta_fai = ref_map["fai"]
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
                ref_dict = ref_map["dict"],
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
                        ref_fasta = ref_map["fasta"],
                        ref_fasta_fai = ref_map["fai"],
                        prefix = prefix,
                        tandem_repeat_bed = tandem_repeat_bed,
                        is_ccs = false,
                        zones = PickGcpZones.zones
                }

                call Utils.InferSampleName {
                    input:
                        bam = SubsetBam.subset_bam,
                        bai = SubsetBam.subset_bai
                }

            }

            call VariantUtils.MergePerChrCalls as MergePBSVVCFs {
                input:
                    vcfs     = RunPBSV.vcf,
                    ref_dict = ref_map["dict"],
                    prefix   = prefix + ".pbsv"
            }

        }

        if (!fast_less_sensitive_sv) {
            call PBSV.RunPBSV as PBSVslow {
                input:
                    bam = bam,
                    bai = bai,
                    ref_fasta = ref_map["fasta"],
                    ref_fasta_fai = ref_map["fai"],
                    prefix = prefix,
                    tandem_repeat_bed = tandem_repeat_bed,
                    is_ccs = false,
                    zones = PickGcpZones.zones
            }

            call VariantUtils.ZipAndIndexVCF as ZipAndIndexPBSV {input: vcf = PBSVslow.vcf }
        }

        call Sniffles2.SampleSV as Sniffles2SV {
            input:
                bam    = bam,
                bai    = bai,
                minsvlen = minsvlen,
                sample_id = sample_id,
                prefix = prefix
        }

        call VariantUtils.ZipAndIndexVCF as ZipAndIndexSnifflesVCF {
            input:
                vcf = Sniffles2SV.vcf
        }
    }

    output {
        File? sniffles_vcf = ZipAndIndexSnifflesVCF.vcfgz
        File? sniffles_tbi = ZipAndIndexSnifflesVCF.tbi

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
