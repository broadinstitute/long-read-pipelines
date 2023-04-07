version 1.0

import "../../../tasks/Utility/FastqUtils.wdl" as FQU

workflow FASTQstats {
    meta {
        desription:
        "Collect some basic stats from a FASTQ file."
    }
    parameter_meta {
        fastq: "file to collect stats on"
        seq_type: "argument to the --seq-type paramter of seqkit"
        exclude_len_threshold: "Sequeces shorter than this will be dropped from analysis."
    }

    input {
        File fastq
        String seq_type = "dna"
        Int? exclude_len_threshold
    }

    if (defined(exclude_len_threshold)) {
        call FQU.FilterByLength   { input: fq = fastq, threshold = select_first([exclude_len_threshold])}
        # call FQU.FilterByLenSeqTk { input: fastq = fastq, exclude_len_threshold = select_first([exclude_len_threshold])}
    }

    call FQU.Stats { input: fastq = select_first([FilterByLength.res, fastq]), seq_type = seq_type }

    output {
        Map[String, Float] stats = Stats.res
    }
}
