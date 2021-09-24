version 1.0

##########################################################################################
# This pipeline calls SVs on an input LR BAM using various known SV algorithms
# that are specifically designed to work with long read data.
# Each individual task/algo. is directly callable, if so desired.
##########################################################################################

import "Structs.wdl"
import "Utils.wdl"
import "VariantUtils.wdl"

import "DeepVariant.wdl" as DV
import "Longshot.wdl"

import "PBSV.wdl"
import "Sniffles.wdl"
import "CuteSV.wdl"
import "SVIM.wdl"

workflow CallVariants {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File? sites_vcf
        File? sites_vcf_tbi
        String? chr

        String prefix

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
    }

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

        call PBSV.Discover {
            input:
                bam               = SubsetBam.subset_bam,
                bai               = SubsetBam.subset_bai,
                ref_fasta         = ref_fasta,
                ref_fasta_fai     = ref_fasta_fai,
                tandem_repeat_bed = tandem_repeat_bed,
                chr               = contig,
                prefix            = prefix
        }

        call PBSV.Call {
            input:
                svsigs        = [ Discover.svsig ],
                ref_fasta     = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                prefix        = prefix
        }



        # Pending a bug fix
#        call DV.PEPPER {
#            input:
#                bam           = SubsetBam.subset_bam,
#                bai           = SubsetBam.subset_bai,
#                ref_fasta     = ref_fasta,
#                ref_fasta_fai = ref_fasta_fai,
#                chr           = contig,
#                preset        = "ONT"
#        }

        call Longshot.Longshot {
            input:
                bam           = SubsetBam.subset_bam,
                bai           = SubsetBam.subset_bai,
                ref_fasta     = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                sites_vcf     = sites_vcf,
                sites_vcf_tbi = sites_vcf_tbi,
                phase         = false,
                chr           = contig
        }
    }
# gather step
    call VariantUtils.MergePerChrCalls as MergePBSVVCFs {
        input:
            vcfs     = Call.vcf,
            ref_dict = ref_dict,
            prefix   = prefix + ".pbsv"
    }


    call VariantUtils.MergePerChrCalls as MergeLongshotVCFs {
        input:
            vcfs     = Longshot.vcf,
            ref_dict = ref_dict,
            prefix   = prefix + ".longshot"
    }



    output {

        File longshot_vcf = MergeLongshotVCFs.vcf
        File longshot_tbi = MergeLongshotVCFs.tbi
        File pbsv_vcf = MergePBSVVCFs.vcf
    }
}