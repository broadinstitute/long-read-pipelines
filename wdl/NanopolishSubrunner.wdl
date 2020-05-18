version 1.0

import "tasks/Nanopolish.wdl" as Nanopolish
import "tasks/Finalize.wdl" as FF

workflow NanopolishSubrunner {
    input {
        Array[File] fast5_and_indexes
        File draft_assembly_fasta
        File draft_assembly_fai
        File draft_alignment_bam
        File draft_alignment_bai

        String output_dir

        File nanopolish_range
    }

    call Nanopolish.NanopolishVariants {
        input:
            fast5_and_indexes = fast5_and_indexes,
            draft_assembly_fasta = draft_assembly_fasta,
            draft_assembly_fai = draft_assembly_fai,
            draft_alignment_bam = draft_alignment_bam,
            draft_alignment_bai = draft_alignment_bai,
            nanopolish_range = nanopolish_range
    }

    call FF.FinalizeToDir {
        input:
            files = NanopolishVariants.vcfs,
            outdir = output_dir
    }
}