version 1.0

##########################################################################################
# This pipeline calls SVs on an input LR BAM using various known SV algorithms
# -- calls the following tasks: Call and Discover in PBSV.wdl
##########################################################################################

import "Structs.wdl"
import "Utils.wdl"
import "VariantUtils.wdl"

import "DeepVariant.wdl" as DV

import "PBSV.wdl"
import "Sniffles.wdl"
import "CuteSV.wdl"
import "SVIM.wdl"

workflow run_pbsv {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

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


    call PBSV.Discover {
        input:
            bam               = bam,
            bai               = bai,
            ref_fasta         = ref_fasta,
            ref_fasta_fai     = ref_fasta_fai,
            tandem_repeat_bed = tandem_repeat_bed,
            prefix            = prefix
    }

    call PBSV.Call {
        input:
            svsigs        = [ Discover.svsig ],
            ref_fasta     = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ccs           = true,
            prefix        = prefix
    }


    output {
        File Discover_out = Discover.svsig
        File Call_out = Call.vcf

    }
}