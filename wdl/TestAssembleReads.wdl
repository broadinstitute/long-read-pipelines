version 1.0

import "Structs.wdl"
import "AssembleReads.wdl" as ASM

workflow TestAssembleReads {
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

    # assemble remaining reads
    call ASM.AssembleReadsWithCanu as AssembleMTFromRemaining {
        input:
            reads = SelectMTFromRemaining.reads,
            target_size = "16.6k",
            platform = platform,
            is_corrected = false,
            prefix = "mt"
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

    # assemble corrected reads
    call ASM.AssembleReadsWithCanu as AssembleMTFromCorrected {
        input:
            reads = SelectMTFromCorrected.reads,
            target_size = "16.6k",
            platform = platform,
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

    # align assembly of remaining reads
    call ASM.AlignContigs as AlignedMTFromRemaining {
        input:
            contigs = AssembleMTFromRemaining.contigs_fasta,
            ref_fasta = ref_fasta,
            SM = "test",
            ID = "remaining",
            PL = "PACBIO",
            is_corrected = true
    }

    # align assembly of corrected reads
    call ASM.AlignContigs as AlignedMTFromCorrected {
        input:
            contigs = AssembleMTFromCorrected.contigs_fasta,
            ref_fasta = ref_fasta,
            SM = "test",
            ID = "corrected",
            PL = "PACBIO",
            is_corrected = true
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

    # call haploid variants on corrected assembly
    call ASM.CallHaploidVariants as CallHaploidVariantsFromCorrected {
        input:
            bam = AlignedMTFromCorrected.aligned_bam,
            bai = AlignedMTFromCorrected.aligned_bai,
            ref_fasta = ref_fasta,
            prefix = "corrected"
    }

    # call haploid variants on remaining assembly
    call ASM.CallHaploidVariants as CallHaploidVariantsFromRemaining {
        input:
            bam = AlignedMTFromRemaining.aligned_bam,
            bai = AlignedMTFromRemaining.aligned_bai,
            ref_fasta = ref_fasta,
            prefix = "remaining"
    }

    # call haploid variants on combined assembly
    call ASM.CallHaploidVariants as CallHaploidVariantsFromCombined {
        input:
            bam = AlignedMTFromCombined.aligned_bam,
            bai = AlignedMTFromCombined.aligned_bai,
            ref_fasta = ref_fasta,
            prefix = "combined"
    }
}

