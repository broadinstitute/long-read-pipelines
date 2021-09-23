version 1.0


##########################################################################################
## A workflow that runs run_pbsv.wdl if given chromosome input and runs
## CallVariantsPBCCS_wholechromosome.wdl if not given chromosome input
##########################################################################################

import "run_pbsv.wdl" as slow

import "CallVariantsPBCCS_wholechromosome.wdl" as fast

workflow chromosome_condition {

    input{
        File bam
        File bai
        File ref_fasta
        File ref_fasta_fai
        File ref_dict
        String prefix
        Boolean fast_less_sensitive
        File? tandem_repeat_bed
        RuntimeAttr? runtime_attr_override
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


    if (!fast_less_sensitive) {
        call slow.run_pbsv as run_pbsv {

        input:

            bam = bam, #assgin input to workflow
            bai = bai,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_dict = ref_dict,
            prefix = prefix,
            tandem_repeat_bed = tandem_repeat_bed

        }

    }

    if (fast_less_sensitive) {

        call fast.CallVariants as run_pbsv2 {

        input:
            bam = bam,
            bai = bai,

            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            tandem_repeat_bed = tandem_repeat_bed,
            ref_dict = ref_dict,
            prefix = prefix
        }


    }

        output {

        File? vcf = if (fast_less_sensitive) then run_pbsv2.pbsv_vcf else run_pbsv.Call_out
        #File svsig = if (fast_less_sensitive) then run_pbsv2.Discover.svsig else run_pbsv.Discover_out

    }

}
