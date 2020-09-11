version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.32/wdl/tasks/SoftClipper.wdl" as SoftClipper
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.32/wdl/tasks/Finalize.wdl" as FF

workflow AssistedSoftClipperRunner {
    input {
        String reads_fastq
        File reference_fasta
        File aid_reference_fasta
        Int rounds
        Int clipping_threshold
        Int ref_conflict_threshold

        File output_dir
    }

    call SoftClipper.SplitSoftClippedReadsAssisted {
        input:
            reads_fastq = reads_fastq,
            reference_fasta = reference_fasta,
            aid_reference_fasta = aid_reference_fasta,
            rounds = rounds,
            clipping_threshold = clipping_threshold,
            ref_conflict_threshold = ref_conflict_threshold,
    }

    call FF.FinalizeToDir {
        input:
            files = flatten([SplitSoftClippedReadsAssisted.split_reads, SplitSoftClippedReadsAssisted.conflicting_alignments]),
            outdir = output_dir
    }
}