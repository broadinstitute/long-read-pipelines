version 1.0

##########################################################################################
# This pipeline calls SVs on an input LR BAM using various known SV algorithms
# that are specifically designed to work with long read data.
# Each individual task/algo. is directly callable, if so desired.
##########################################################################################

import "Structs.wdl"
import "PBSV.wdl"
import "Utils.wdl"
import "VariantUtils.wdl"
import "Sniffles.wdl"

workflow CallVariants {
    input {
        File bam
        File bai
        File ref_fasta
        File ref_fasta_fai
        File ref_dict
        String prefix
        File? tandem_repeat_bed
        Boolean fast_less_sensitive
    }

    parameter_meta {
        bam:               "input BAM from which to call SVs"
        bai:               "index accompanying the BAM"
        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"
        ref_dict:          "sequence dictionary accompanying the reference"
        prefix:            "prefix for output files"
        tandem_repeat_bed: "BED file containing TRF finder (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"
    }

    if (fast_less_sensitive) {

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
                    bam = bam,
                    bai = bai,
                    ref_fasta = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,
                    ref_dict = ref_dict,
                    prefix = prefix,
                    tandem_repeat_bed = tandem_repeat_bed,
                    ccs = true
            }

            call Sniffles.Sniffles {
                input:
                    bam    = SubsetBam.subset_bam,
                    bai    = SubsetBam.subset_bai,
                    chr    = contig,
                    prefix = prefix
            }

    }

        call VariantUtils.MergePerChrCalls as MergePBSVVCFs {
        input:
            vcfs     = RunPBSV.vcf,
            ref_dict = ref_dict,
            prefix   = prefix + ".pbsv"
        }

        call VariantUtils.MergePerChrCalls as MergeSnifflesVCFs {
        input:
            vcfs     = Sniffles.vcf,
            ref_dict = ref_dict,
            prefix   = prefix + ".sniffles"
        }

    }

    if (!fast_less_sensitive) {

        call PBSV.RunPBSV as PBSVslow {
            input:
                bam = bam,
                bai = bai,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_dict = ref_dict,
                prefix = prefix,
                tandem_repeat_bed = tandem_repeat_bed,
                ccs = true
        }


        call Sniffles.Sniffles as SnifflesSlow {
            input:
                bam    = bam,
                bai    = bai,
                prefix = prefix
        }

    }


    output {
        File sniffles_vcf = select_first([MergeSnifflesVCFs.vcf, SnifflesSlow.vcf])
        File pbsv_vcf = select_first([MergePBSVVCFs.vcf, PBSVslow.vcf])
    }

}


