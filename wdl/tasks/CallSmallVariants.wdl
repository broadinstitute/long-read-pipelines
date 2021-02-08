version 1.0

##########################################################################################
# This pipeline calls small variants (currently only SNVs) on an input LR BAM using
# known algorithms that are specifically designed to work with long read data.
##########################################################################################

import "Structs.wdl"
import "Utils.wdl" as Utils
import "Longshot.wdl" as Longshot
import "DeepVariant.wdl" as DV

workflow CallSmallVariants {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String preset
    }

    parameter_meta {
        bam: "input BAM from which to call SNVs"
        bai: "index accompanying the BAM"

        ref_fasta:     "reference to which the BAM was aligned to"
        ref_fasta_fai: "index accompanying the reference"
        ref_dict:      "dictionary accompanying the reference"
    }

    call Utils.MakeChrIntervalList { input: ref_dict = ref_dict }

    scatter (chr_info in MakeChrIntervalList.chrs) {
        call Longshot.Longshot {
            input:
                bam = bam,
                bai = bai,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                chr = chr_info[0]
        }

        call DV.DeepVariant {
            input:
                bam           = bam,
                bai           = bai,
                ref_fasta     = ref_fasta,
                ref_fai       = ref_fasta_fai,
                model_class   = "PACBIO",
                chr           = chr_info[0]
        }
    }

    call Longshot.MergeLongshotCalls {
        input:
            vcfs = Longshot.vcf,
            ref_dict = ref_dict,
            prefix = basename(bam, ".bam")
    }

    output {
        File longshot_vcf = MergeLongshotCalls.vcf
        File longshot_tbi = MergeLongshotCalls.tbi
#        File vcf = DeepVariant.vcf
#        File vcf_tbi = DeepVariant.vcf_tbi
#        File gvcf = DeepVariant.gvcf
#        File gvcf_tbi = DeepVariant.gvcf_tbi
    }
}
