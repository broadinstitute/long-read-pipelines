version 1.0

##########################################################################################
# This pipeline calls SVs on an input LR BAM using various known SV algorithms
# that are specifically designed to work with long read data.
# Each individual task/algo. is directly callable, if so desired.
##########################################################################################

import "Utils.wdl"
import "VariantUtils.wdl"

import "PBSV.wdl"
import "Sniffles2.wdl" as Sniffles2


workflow CallVariants {
    meta {
        descrition: "A workflow for calling small and/or structural variants from an aligned CLR BAM file."
    }
    input {
        File bam
        File bai
        Int minsvlen = 50
        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String prefix
        String sample_id

        Boolean call_svs
        Boolean pbsv_call_per_chr
        File? tandem_repeat_bed
    }

    parameter_meta {
        bam:               "input BAM from which to call SVs"
        bai:               "index accompanying the BAM"

        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"
        ref_dict:          "sequence dictionary accompanying the reference"

        prefix:            "prefix for output files"

        tandem_repeat_bed: "BED file containing TRF finder (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"
        minsvlen:       "Minimum SV length in bp (default: 50)"
        pbsv_call_per_chr:      "when using PBSV, make calls per chromosome then merge (trade lower sensitivity for faster speed)"
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
        if (pbsv_call_per_chr) {

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
            }

            call VariantUtils.MergePerChrCalls as MergePBSVVCFs {
                input:
                    vcfs     = RunPBSV.vcf,
                    ref_dict = ref_dict,
                    prefix   = prefix + ".pbsv"
            }
        }

        if (!pbsv_call_per_chr) {

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
        File? sniffles_vcf = select_first([ZipAndIndexSnifflesVCF.vcfgz])
        File? sniffles_tbi = select_first([ZipAndIndexSnifflesVCF.tbi])
        File? sniffles_snf = select_first([Sniffles2SV.snf])

        File? pbsv_vcf = select_first([MergePBSVVCFs.vcf, ZipAndIndexPBSV.vcfgz])
        File? pbsv_tbi = select_first([MergePBSVVCFs.tbi, ZipAndIndexPBSV.tbi])
    }
}
