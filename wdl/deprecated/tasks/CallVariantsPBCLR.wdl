version 1.0

##########################################################################################
# This pipeline calls SVs on an input LR BAM using various known SV algorithms
# that are specifically designed to work with long read data.
# Each individual task/algo. is directly callable, if so desired.
##########################################################################################

import "../../tasks/Utility/Utils.wdl"
import "../../tasks/Utility/VariantUtils.wdl"

import "../../tasks/VariantCalling/PBSV.wdl"
import "Sniffles.wdl"


workflow CallVariants {
    meta {
        descrition: "A workflow for calling small and/or structural variants from an aligned CLR BAM file."
    }

    parameter_meta {
        bam:               "input BAM from which to call SVs"
        bai:               "index accompanying the BAM"

        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"
        ref_dict:          "sequence dictionary accompanying the reference"

        call_small_variants: "if true, will attempt to call small variants"
        call_small_vars_on_mitochondria: "if false, will not attempt to call variants on mitochondria"
        fast_less_sensitive_sv: "if true, will run PBSV in a less sensitive mode, which is faster"
        sites_vcf:     "for use with Clair, the small variant caller"
        sites_vcf_tbi: "for use with Clair, the small variant caller"

        prefix:            "prefix for output files"

        tandem_repeat_bed: "BED file containing TRF finder (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"
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

        Boolean call_small_variants = false
        Boolean call_small_vars_on_mitochondria
    }

    call Utils.RandomZoneSpewer as arbitrary {input: num_of_zones = 3}

    ######################################################################
    # Block for small variants handling
    ######################################################################
    # todo: use NanoCaller, Clair isn't going to support CLR
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

                call PBSV.Discover as pbsv_discover_chr {
                    input:
                        bam = SubsetBam.subset_bam,
                        bai = SubsetBam.subset_bai,
                        is_hifi = true,
                        ref_fasta = ref_fasta,
                        ref_fasta_fai = ref_fasta_fai,
                        tandem_repeat_bed = tandem_repeat_bed,
                        chr = contig_for_sv,
                        prefix = prefix,
                        zones = arbitrary.zones
                }

                call Sniffles.Sniffles {
                    input:
                        bam    = SubsetBam.subset_bam,
                        bai    = SubsetBam.subset_bai,
                        chr    = contig_for_sv,
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
            call PBSV.Call as pbsv_wg_call {
                input:
                    svsigs = pbsv_discover_chr.svsig,
                    ref_fasta = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,
                    is_hifi = true,
                    prefix = prefix + ".pbsv",
                    zones = arbitrary.zones
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
                    is_hifi = true,
                    zones = arbitrary.zones
            }

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

        File? pbsv_vcf = select_first([pbsv_wg_call.vcf, PBSVslow.vcf])
        File? pbsv_tbi = select_first([pbsv_wg_call.tbi, PBSVslow.tbi])
    }
}
