version 1.0

import "tasks/SoftClipper.wdl" as SoftClipper
import "tasks/Finalize.wdl" as FF

workflow SoftClipperAidedRunner {
    input {
        String reads_fastq
        File reference_fasta
        File aid_reference_fasta
        Int rounds
        Int clipping_threshold
        Boolean keep_all_ref_conflicts

        File output_dir
    }

    call SoftClipper.SplitSoftClippedReadsAided {
        input:
            reads_fastq = reads_fastq,
            reference_fasta = reference_fasta,
            aid_reference_fasta = aid_reference_fasta,
            rounds = rounds,
            clipping_threshold = clipping_threshold,
            keep_all_ref_conflicts = keep_all_ref_conflicts
    }

    call FF.FinalizeToDir {
        input:
            files = SplitSoftClippedReadsAided.split_reads,
            outdir = output_dir
    }
}