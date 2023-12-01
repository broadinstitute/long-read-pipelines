version 1.0

import "../../../tasks/Utility/FastqUtils.wdl" as FQU
import "../../../tasks/Utility/BAMutils.wdl" as BU

workflow FASTQstats {
    meta {
        desription:
        "Collect some basic stats from a FASTQ (or BAM) file."
    }
    parameter_meta {
        reads: "file to collect stats on"
        file_type: "type of file: accepted values are [FASTQ, BAM] (regardless of gz, bgz)"
        seq_type: "argument to the --seq-type paramter of seqkit"
        exclude_len_threshold: "Sequeces shorter than this will be dropped from analysis; no effect if not provided"
    }

    input {
        File reads
        String file_type
        String seq_type = "dna"
        Int? exclude_len_threshold
    }

    output {
        Map[String, Float] stats = Stats.res
    }

    if ('BAM' == file_type) {
        call BU.BamToFastq { input: bam = reads, prefix = basename(reads, ".bam"), disk_type = 'SSD' }
    }
    File formatted_input = select_first([BamToFastq.reads_fq, reads])
    if (defined(exclude_len_threshold)) {
        call FQU.FilterByLength   { input: fq = formatted_input, threshold = select_first([exclude_len_threshold])}
        # call FQU.FilterByLenSeqTk { input: fastq = fastq, exclude_len_threshold = select_first([exclude_len_threshold])}
    }
    File filtered_input = select_first([FilterByLength.res, formatted_input])

    call FQU.Stats { input: fastq = filtered_input, seq_type = seq_type }
}
