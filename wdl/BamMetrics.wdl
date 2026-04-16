version 1.0

workflow BamMetrics {
    meta {
        description: "Minimal workflow to compute coverage and read-length metrics from an aligned BAM."
    }

    input {
        File bam
        File bai
        File ref_fasta
    }

    call ComputeGenomeLength {
        input:
            fasta = ref_fasta
    }

    call NanoPlotFromBam {
        input:
            bam = bam,
            bai = bai
    }

    output {
        Float aligned_coverage = NanoPlotFromBam.stats_map["number_of_bases_aligned"] / ComputeGenomeLength.length

        Float read_length_mean = NanoPlotFromBam.stats_map["mean_read_length"]
        Float read_length_median = NanoPlotFromBam.stats_map["median_read_length"]
        Float read_length_stdev = NanoPlotFromBam.stats_map["read_length_stdev"]
        Float read_length_n50 = NanoPlotFromBam.stats_map["n50"]

        Float aligned_num_reads = NanoPlotFromBam.stats_map["number_of_reads"]
        Float aligned_num_bases = NanoPlotFromBam.stats_map["number_of_bases_aligned"]
        Float aligned_fraction_bases = NanoPlotFromBam.stats_map["fraction_bases_aligned"]

        Float average_identity = NanoPlotFromBam.stats_map["average_identity"]
        Float median_identity = NanoPlotFromBam.stats_map["median_identity"]

        File nanoplot_stats = NanoPlotFromBam.stats
    }
}

task ComputeGenomeLength {
    input {
        File fasta
    }

    Int disk_gb = 2 * ceil(size(fasta, "GB"))

    command <<<
        set -euxo pipefail

        samtools dict ~{fasta} \
            | grep '^@SQ' \
            | awk '{ print $3 }' \
            | sed 's/LN://' \
            | awk '{ sum += $1 } END { print sum }' \
            > genome_length.txt
    >>>

    output {
        Float length = read_float("genome_length.txt")
    }

    runtime {
        cpu: 1
        memory: "1 GiB"
        disks: "local-disk " + disk_gb + " HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
    }
}

task NanoPlotFromBam {
    input {
        File bam
        File bai
    }

    Int disk_gb = 2 * ceil(size(bam, "GB")) + 10

    command <<<
        set -euxo pipefail

        touch ~{bai}

        num_core=$(grep -c '^processor' /proc/cpuinfo)

        NanoPlot -t ${num_core} \
                 -c orangered \
                 --N50 \
                 --tsv_stats \
                 --no_supplementary \
                 --verbose \
                 --bam "~{bam}"

        cat NanoStats.txt \
            | grep -v -e '^Metrics' -e '^highest' -e '^longest' \
            | sed 's/ >/_/' \
            | sed 's/://' \
            | awk '{ print $1 "\t" $2 }' \
            > map.txt
    >>>

    output {
        File stats = "NanoStats.txt"
        Map[String, Float] stats_map = read_map("map.txt")
    }

    runtime {
        cpu: 8
        memory: "24 GiB"
        disks: "local-disk " + disk_gb + " LOCAL"
        docker: "us.gcr.io/broad-dsp-lrma/lr-nanoplot:1.40.0-1"
    }
}
