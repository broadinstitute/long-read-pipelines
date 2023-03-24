version 1.0

##########################################################################################
# This pipeline calls SVs on an input LR BAM using various known SV algorithms
# that are specifically designed to work with long read data.
# Each individual task/algo. is directly callable, if so desired.
##########################################################################################

import "../Utility/Utils.wdl"
import "../Utility/VariantUtils.wdl"

import "PBSV.wdl"
import "Sniffles.wdl"


workflow CallVariants {
    meta {
        descrition: "A workflow for calling small and/or structural variants from an aligned CLR BAM file."
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

    parameter_meta {
        bam:               "input BAM from which to call SVs"
        bai:               "index accompanying the BAM"

        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"
        ref_dict:          "sequence dictionary accompanying the reference"

        call_small_vars_on_mitochondria: "if false, will not attempt to call variants on mitochondria"
        sites_vcf:     "for use with Clair, the small variant caller"
        sites_vcf_tbi: "for use with Clair, the small variant caller"

        prefix:            "prefix for output files"

        tandem_repeat_bed: "BED file containing TRF finder (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"
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
    }
}
