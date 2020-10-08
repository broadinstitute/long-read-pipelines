version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.35/wdl/tasks/Canu.wdl" as Canu
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.35/wdl/tasks/SoftClipper.wdl" as SoftClipper
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.35/wdl/tasks/Finalize.wdl" as FF

workflow SplitSoftClipsAndAssemble {
    input {
        String reads_fastq
        File reference_fasta
        Int rounds_of_splitting
        Int clipping_threshold

        String canu_output_file_prefix
        String canu_genome_size
        Float canu_correctedErrorRate

        File out_dir
    }

    call SoftClipper.SplitSoftClippedReads {
        input:
            reads_fastq = reads_fastq,
            reference_fasta = reference_fasta,
            rounds = rounds_of_splitting,
            clipping_threshold = clipping_threshold
    }

    call Canu.CorrectTrimAssemble {
        input:
            output_file_prefix = canu_output_file_prefix,
            genome_size = canu_genome_size,
            reads = SplitSoftClippedReads.most_split_read,
            correct_corrected_error_rate = canu_correctedErrorRate,
            trim_corrected_error_rate = canu_correctedErrorRate,
            assemble_corrected_error_rate = canu_correctedErrorRate,
    }

    call FF.FinalizeToDir as FinalizeAssembly {
        input:
            files = [CorrectTrimAssemble.canu_contigs_fasta],
            outdir = out_dir
    }
}