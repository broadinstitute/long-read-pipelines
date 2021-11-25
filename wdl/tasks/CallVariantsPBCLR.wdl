version 1.0

##########################################################################################
# This pipeline calls SVs on an input LR BAM using various known SV algorithms
# that are specifically designed to work with long read data.
# Each individual task/algo. is directly callable, if so desired.
##########################################################################################

import "Structs.wdl"
import "Utils.wdl"
import "VariantUtils.wdl"

import "PBSV.wdl"
import "Sniffles.wdl"

import "Clair.wdl" as Clair3

workflow CallVariants {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String prefix

        Boolean fast_less_sensitive_sv

        Boolean call_small_vars_on_mitochondria
        File? sites_vcf
        File? sites_vcf_tbi

        File? tandem_repeat_bed
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

    ######################################################################
    # Block for small variants handling
    ######################################################################
    Array[String] default_filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']
    Array[String] use_filter = if (call_small_vars_on_mitochondria) then default_filter else flatten([['chrM'],default_filter])
    call Utils.MakeChrIntervalList as SmallVariantsScatterPrepp {
        input:
            ref_dict = ref_dict,
            filter = use_filter
    }

    scatter (c in SmallVariantsScatterPrepp.chrs) {
        String contig = c[0]

        call Utils.SubsetBam as SmallVariantsScatter {
            input:
                bam = bam,
                bai = bai,
                locus = contig
        }

        call Clair3.Clair {
            input:
                bam = SmallVariantsScatter.subset_bam,
                bai = SmallVariantsScatter.subset_bai,

                ref_fasta     = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,

                sites_vcf = sites_vcf,
                sites_vcf_tbi = sites_vcf_tbi,

                preset = "CLR"
        }
    }

    call VariantUtils.MergeAndSortVCFs as MergeAndSortClairVCFs {
        input:
            vcfs = Clair.vcf,
            tbis = Clair.vcf_tbi,
            ref_fasta_fai = ref_fasta_fai,
            prefix = prefix + ".clair"
    }

    call VariantUtils.MergeAndSortVCFs as MergeAndSortClair_gVCFs {
        input:
            vcfs = Clair.gvcf,
            tbis = Clair.gvcf_tbi,
            ref_fasta_fai = ref_fasta_fai,
            prefix = prefix + ".clair.g"
    }

    ######################################################################
    # Block for SV handling
    ######################################################################
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
                    is_ccs = false
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

        call VariantUtils.MergeAndSortVCFs as MergeSnifflesVCFs {
            input:
                vcfs   = FixSnifflesVCF.sortedVCF,
                tbis   = FixSnifflesVCF.tbi,
                ref_fasta_fai = ref_fasta_fai,
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
                is_ccs = false
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

    output {
        File sniffles_vcf = select_first([MergeSnifflesVCFs.vcf, ZipAndIndexSniffles.sortedVCF])
        File sniffles_tbi = select_first([MergeSnifflesVCFs.tbi, ZipAndIndexSniffles.tbi])

        File pbsv_vcf = select_first([MergePBSVVCFs.vcf, ZipAndIndexPBSV.vcfgz])
        File pbsv_tbi = select_first([MergePBSVVCFs.tbi, ZipAndIndexPBSV.tbi])

        File clair_vcf = MergeAndSortClairVCFs.vcf
        File clair_tbi = MergeAndSortClairVCFs.tbi

        File clair_gvcf = MergeAndSortClair_gVCFs.vcf
        File clair_gtbi = MergeAndSortClair_gVCFs.tbi
    }
}
