version 1.0

##########################################################################################
# A workflow that runs Translocator
##########################################################################################

import "tasks/Translocator.wdl" as Tr

workflow Translocator {
    input {
        File aligned_bam
        File ref_fasta

        String prefix
        Float? min_het_af
    }
    call Tr.Translocator {input: aligned_bam = aligned_bam, ref_fasta = ref_fasta, prefix = prefix, min_het_af_override = min_het_af}
}
