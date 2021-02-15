version 1.0

##########################################################################################
# This pipeline calls SVs on an input LR BAM using various known SV algorithms
# that are specifically designed to work with long read data.
# Each individual task/algo. is directly callable, if so desired.
##########################################################################################

import "Structs.wdl"
import "CuteSV.wdl"
import "Sniffles.wdl"
import "PBSV.wdl"
import "SVIM.wdl"

workflow CallSVsONT {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai
        File? tandem_repeat_bed
    }

    parameter_meta {
        bam: "input BAM from which to call SVs"
        bai: "index accompanying the BAM"

        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"
        tandem_repeat_bed: "BED file containing TRF finder (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"
    }

#    call PBSV.PBSV {
#        input:
#            bam = bam,
#            bai = bai,
#            ref_fasta = ref_fasta,
#            ref_fasta_fai = ref_fasta_fai,
#            tandem_repeat_bed = tandem_repeat_bed,
#            prefix = basename(bam, ".bam")
#    }

    call Sniffles.Sniffles {
        input:
            bam = bam,
            bai = bai,
            prefix = basename(bam, ".bam")
    }

    call SVIM.SVIM {
        input:
            bam = bam,
            bai = bai,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            prefix = basename(bam, ".bam")
    }

#    call CuteSV.CuteSV {
#        input:
#            bam = bam,
#            bai = bai,
#            ref_fasta = ref_fasta,
#            prefix = basename(bam, ".bam"),
#            preset = "ont"
#    }

    output {
#        File pbsv_vcf = PBSV.vcf
        File sniffles_vcf = Sniffles.vcf
        File svim_vcf = SVIM.vcf
#        File cutesv_vcf = CuteSV.vcf
    }
}
