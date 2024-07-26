version 1.0

workflow MosdepthCoverageStats {

    meta {
        description: "Calculate coverage statistics using mosdepth."
    }
    parameter_meta {
        aligned_bam: "Aligned BAM file."
        aligned_bai: "Aligned BAM index file."
        bed_file: "BED file containing regions of interest."
        bin_length: "Length of bins to use for coverage calculation."
        preemptible_tries: "Number of times to retry a preempted task."
    }

    input {
        File aligned_bam
        File aligned_bai
        File? bed_file

        Int bin_length = 1000

        # Runtime parameters
        Int preemptible_tries = 3
    }

    call MosDepthOverBed {
        input:
            bam = aligned_bam,
            bai = aligned_bai,
            bed = bed_file,
            bin_length = bin_length,
            preemptible_tries = preemptible_tries,
    }

    String basename = basename(aligned_bam, ".bam")

    call CoverageStats {
        input:
            mosdepth_regions = MosDepthOverBed.regions,
            basename_input = basename,
            preemptible_tries = preemptible_tries,
    }

    output {
        File cov_stat_summary_file      = CoverageStats.cov_stat_summary_file
        Map[String, Float] cov_stat_summary = CoverageStats.cov_stat_summary
    }

}


task MosDepthOverBed {

    meta {
        description: "Calculate coverage using mosdepth."
    }
    parameter_meta {
        bam: "Aligned BAM file."
        bai: "Aligned BAM index file."
        bed: "BED file containing regions of interest."
        bin_length: "Length of bins to use for coverage calculation. default is 1000. If bed file is provided, this parameter is ignored."
        chrom: "Chromosome to calculate coverage for."
        no_per_base: "Do not calculate per-base coverage."
        fast_mode: "Use fast mode."
        mapq: "Minimum mapping quality to consider."
        thresholds: "Comma-separated list of thresholds to use for coverage calculation. (e.g. 1,10,20,30)."
        preemptible_tries: "Number of times to retry a preempted task."
    }

    input {
        File bam
        File bai
        File? bed
        Int threads = 4
        String? chrom
        Int? bin_length
        String? thresholds
        Boolean no_per_base = true
        Boolean fast_mode = true
        Int mapq = 1

        # Runtime parameters
        Int preemptible_tries = 3
    }

    Int disk_size = 2 * ceil(size(bam, "GB") + size(bai, "GB"))
    String basename = basename(bam, ".bam")
    String prefix = "~{basename}.coverage_over_bed"

    command {
        set -euxo pipefail

        # Create symbolic links for bam and bai in the current working directory
        ln -s ~{bam} ./~{basename}.bam
        ln -s ~{bai} ./~{basename}.bai

        mosdepth \
        ~{true="-n" false="" no_per_base} \
        ~{true="-x" false="" fast_mode} \
        -t ~{threads} \
        ~{"-c " + chrom} \
        ~{"-b " + select_first([bed, bin_length])} \
        ~{"-Q " + mapq} \
        ~{"-T " + thresholds} \
        ~{prefix} ./~{basename}.bam
    }

    output {
        File global_dist      = "~{prefix}.mosdepth.global.dist.txt"
        File region_dist      = "~{prefix}.mosdepth.region.dist.txt"
        File regions          = "~{prefix}.regions.bed.gz"
        File regions_csi      = "~{prefix}.regions.bed.gz.csi"
    }

    runtime {
        cpu:                    4
        memory:                 8 + " GiB"
        disks: "local-disk " +  disk_size + " HDD"
        preemptible:            preemptible_tries
        maxRetries:             1
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-mosdepth:latest"
    }
}

task CoverageStats {

    meta {
        description: "Calculate coverage statistics from mosdepth output."
    }
    parameter_meta {
        mosdepth_regions: "Mosdepth output file containing coverage values."
        cov_col: "Column holding the coverage values."
        basename_input: "Basename to use for output files."
        preemptible_tries: "Number of times to retry a preempted task."
    }

    input {
        File mosdepth_regions
        Int cov_col = 4 # column holding the coverage values
        String? basename_input

        # Runtime parameters
        Int preemptible_tries = 3
    }
    Int round = 2
    String header_suffix = "_coverage"

    Int disk_size = 2*ceil(size(mosdepth_regions, "GB"))
    String basename = select_first([basename_input, basename(mosdepth_regions)])
    String prefix = "~{basename}.coverage_over_bed"

    command {
        set -euxo pipefail

        python3 coverage_stats.py \
        --cov_col ~{cov_col} \
        --round ~{round} \
        --output_prefix ~{prefix} \
        ~{mosdepth_regions}
    }

    output {
        File cov_stat_summary_file = "~{prefix}.cov_stat_summary.json"
        Map[String, Float] cov_stat_summary = read_json("~{prefix}.cov_stat_summary.json")
    }

    runtime {
        cpu:                    2
        memory:                 4 + " GiB"
        disks: "local-disk " +  disk_size + " HDD"
        preemptible:            preemptible_tries
        maxRetries:             1
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-mosdepth:latest"
    }
}
