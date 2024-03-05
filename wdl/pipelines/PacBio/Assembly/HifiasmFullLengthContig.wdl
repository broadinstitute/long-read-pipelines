version 1.0

import "../../../tasks/Assembly/Utility.wdl" as U

workflow HifiasmFullLengthContig{
    input{
        File wholegenomebam
        File wholegenomebai
        File reference_fasta
        String asmregion
        String prefix
        Int nthreads
    }
    call U.extract_reads{input: bam_input=wholegenomebam, bam_index=wholegenomebai, region=asmregion, pref=prefix}
    call U.hifiasm_asm{input: reads=extract_reads.local_fq, prefix=prefix}
    call U.realign_reads{input: wholebam=wholegenomebam, wholebai=wholegenomebai ,assembly=hifiasm_asm.assembly_primary, num_threads=nthreads, pref=prefix}
    call U.mergefq{input: first_round_reads=extract_reads.local_fq, reads=realign_reads.realignedreadsfq, prefix=prefix}
    call U.hifiasm_asm as second_round{input: reads = mergefq.merged_fq, prefix = prefix}


    call U.ragtag_construct as ragtag_hap1{input: asm = second_round.assembly_hap1, ref = reference_fasta, outputfolder = prefix}
    call U.ragtag_construct as ragtag_hap2{input: asm = second_round.assembly_hap2, ref = reference_fasta, outputfolder = prefix}
    call U.extract_bam{input: bam_input=wholegenomebam, bam_index=wholegenomebai, region= asmregion, pref=prefix}
    call U.tsggapcloser_gapfilling as tgs_hap1{input: scaffold_input = ragtag_hap1.scaffold, read_fq = extract_bam.local_fa, outputprefix = prefix}
    call U.tsggapcloser_gapfilling as tgs_hap2{input: scaffold_input = ragtag_hap2.scaffold, read_fq = extract_bam.local_fa, outputprefix = prefix}
    

    meta{
        Purpose:"Local assembly using hifiasm"
    }
    output{
        File ragtag_scaffold_hap1 = ragtag_hap1.scaffold
        File ragtag_scaffold_hap2 = ragtag_hap2.scaffold
        File tgs_scaffold_hap1 = tgs_hap1.scaffold
        File tgs_scaffold_hap2 = tgs_hap2.scaffold
    }
}

