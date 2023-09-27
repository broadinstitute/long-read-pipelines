version 1.0

import "../../../tasks/Utility/FastqUtils.wdl" as FQU
import "../../../tasks/Utility/BAMutils.wdl" as BU

workflow FASTQstatsFromBam {
    meta {
        desription:
        "Collect some basic stats from a BAM file."
    }
    parameter_meta {
        bam: "file to collect stats on"
        seq_type: "argument to the --seq-type paramter of seqkit"
    }

    input {
        File  bam
        File? bai
        String seq_type = "dna"
        String disk_type = "HDD"
    }

    call BU.BamToFastq { input: bam = bam, prefix = basename(bam, ".bam"), disk_type = disk_type }
    call FQU.Stats { input: fastq = BamToFastq.reads_fq, seq_type = seq_type }

    output {
        Map[String, Float] stats = Stats.res
    }
}
