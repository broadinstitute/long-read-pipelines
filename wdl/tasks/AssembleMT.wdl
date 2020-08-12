version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/Structs.wdl"
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/AssembleTarget.wdl" as AT

workflow AssembleMT {
    input {
        File aligned_bam
        File aligned_bai

        String platform

        File ref_fasta
        File ref_fasta_fai
        File ref_dict
        String mt_chr_name

        String prefix
    }

    # select reads from MT
    call AT.SelectReadsFromRegion as SelectMT {
        input:
            bam = aligned_bam,
            bai = aligned_bai,
            region = [ mt_chr_name ]
    }

    # assemble reads
    call AT.AssembleReadsWithCanu as AssembleMT {
        input:
            reads = SelectMT.reads,
            target_size = "16.6k",
            platform = platform,
            is_corrected = if platform == "PACBIO" then true else false,
            prefix = prefix
    }

    # align assembly of combined reads
    call AT.AlignContigs as AlignMTAssembly {
        input:
            contigs = AssembleMT.contigs_fasta,
            ref_fasta = ref_fasta,
            SM = "test",
            ID = "combined",
            PL = platform,
            is_corrected = if platform == "PACBIO" then true else false,
            prefix = prefix
    }

    # call haploid variants on combined assembly
    call AT.CallHaploidVariants as CallHaploidVariants {
        input:
            bam = AlignMTAssembly.aligned_bam,
            bai = AlignMTAssembly.aligned_bai,
            ref_fasta = ref_fasta,
            prefix = prefix
    }

    output {
        File report          = AssembleMT.report

        File contigs_fasta   = AssembleMT.contigs_fasta
        File unassembled     = AssembleMT.unassembled
        File unitigs_fasta   = AssembleMT.unitigs_fasta

        File contigs_layout  = AssembleMT.contigs_layout
        File unitigs_layout  = AssembleMT.unitigs_layout
        File unitigs_bed     = AssembleMT.unitigs_bed

        File contigs_gfa     = AssembleMT.contigs_gfa
        File unitigs_gfa     = AssembleMT.unitigs_gfa

        File mt_aligned_bam  = AlignMTAssembly.aligned_bam
        File mt_aligned_bai  = AlignMTAssembly.aligned_bai

        File mt_calls        = CallHaploidVariants.calls
    }
}

