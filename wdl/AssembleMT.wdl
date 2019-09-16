version 1.0

import "Structs.wdl"
import "AssembleReads.wdl" as ASM

workflow AssembleMT {
    input {
        File corrected_bam
        File corrected_bai
        File remaining_bam
        File remaining_bai

        String platform

        File ref_fasta
        File ref_fasta_fai
        File ref_dict
    }

    # select corrected reads from MT
    call ASM.SelectReadsFromRegion as SelectMTFromCorrected {
        input:
            bam = corrected_bam,
            bai = corrected_bai,
            region = "chrM"
    }

    # select remaining reads from MT
    call ASM.SelectReadsFromRegion as SelectMTFromRemaining {
        input:
            bam = remaining_bam,
            bai = remaining_bai,
            region = "chrM"
    }

    # correct/trim remaining reads
    call ASM.CorrectAndTrimReadsWithCanu as CorrectMTFromRemaining {
        input:
            reads = SelectMTFromRemaining.reads,
            target_size = "16.6k",
            platform = platform,
            is_corrected = false,
            prefix = "mt"
    }

    # combine corrected and remaining reads
    call ASM.CombineReads as CombineReads {
        input:
            reads = [ SelectMTFromCorrected.reads, CorrectMTFromRemaining.trimmed_reads ]
    }

    # assemble combined reads
    call ASM.AssembleReadsWithCanu as AssembleMTFromCombined {
        input:
            reads = CombineReads.reads,
            target_size = "16.6k",
            platform = platform,
            prefix = "mt"
    }

    # align assembly of combined reads
    call ASM.AlignContigs as AlignedMTFromCombined {
        input:
            contigs = AssembleMTFromCombined.contigs_fasta,
            ref_fasta = ref_fasta,
            SM = "test",
            ID = "combined",
            PL = "PACBIO",
            is_corrected = true
    }

    # call haploid variants on combined assembly
    call ASM.CallHaploidVariants as CallHaploidVariantsFromCombined {
        input:
            bam = AlignedMTFromCombined.aligned_bam,
            bai = AlignedMTFromCombined.aligned_bai,
            ref_fasta = ref_fasta,
            prefix = "combined"
    }

    output {
        File report          = AssembleMTFromCombined.report

        File contigs_fasta   = AssembleMTFromCombined.contigs_fasta
        File unassembled     = AssembleMTFromCombined.unassembled
        File unitigs_fasta   = AssembleMTFromCombined.unitigs_fasta

        File contigs_layout  = AssembleMTFromCombined.contigs_layout
        File unitigs_layout  = AssembleMTFromCombined.unitigs_layout
        File unitigs_bed     = AssembleMTFromCombined.unitigs_bed

        File contigs_gfa     = AssembleMTFromCombined.contigs_gfa
        File unitigs_gfa     = AssembleMTFromCombined.unitigs_gfa

        File aligned_bam     = AlignedMTFromCombined.aligned_bam
        File aligned_bai     = AlignedMTFromCombined.aligned_bai
        File calls           = CallHaploidVariantsFromCombined.calls
    }
}

