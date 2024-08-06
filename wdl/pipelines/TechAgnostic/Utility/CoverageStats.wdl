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
        stats_per_interval: "If true, calculate coverage statistics per interval. Must provide a bed file."
    }

    input {
        File aligned_bam
        File aligned_bai
        File? bed_file

        Int? bin_length

        Boolean stats_per_interval = false

        # Runtime parameters
        Int preemptible_tries = 3
        Int mosdepth_mem = 8
    }

    if (stats_per_interval ) {
        if (!defined(bed_file)) {
            call FailWorkflow {
                input:
                    message = "stats_per_interval is set to true, but no bed file is provided."
            }
        }
        scatter ( bed_line in read_lines(select_first([bed_file])) ) {
            call MosDepthOverBed as MosDepthPerInterval {
                input:
                    bam = aligned_bam,
                    bai = aligned_bai,
                    bed = write_lines([bed_line]),
                    preemptible_tries = preemptible_tries,
                    mem = mosdepth_mem,
            }
        }
        Array[File] cov_stat_summary_files  = MosDepthPerInterval.cov_stat_summary_file
        Array[Map[String, Float]] cov_stat_summaries = MosDepthPerInterval.cov_stat_summary
    }



    if (defined(bed_file) && defined(bin_length)) {
        call BinBed {
            input:
                bed_file = bed_file,
                bin_length = bin_length
        }
    }

    # Runs in cases where bed file is provided or bed and bin_length are provided
    if (defined(bed_file) || (defined(bin_length) && defined(BinBed.binned_bed))) {
        call MosDepthOverBed {
            input:
                bam = aligned_bam,
                bai = aligned_bai,
                bed = select_first([BinBed.binned_bed, bed_file]),
                bin_length = bin_length,
                preemptible_tries = preemptible_tries,
                mem = mosdepth_mem,
        }
    }

    # Runs in cases where no bed file is provided, will run regardless of bin_length
    if (!defined(bed_file) ) {
        call MosDepthOverBed as MosDepthNoBed {
            input:
                bam = aligned_bam,
                bai = aligned_bai,
                bin_length = bin_length,
                preemptible_tries = preemptible_tries,
                mem = mosdepth_mem,
        }
    }


    output {
        File cov_stat_summary_file      = select_first([MosDepthOverBed.cov_stat_summary_file, MosDepthNoBed.cov_stat_summary_file])
        Map[String, Float] cov_stat_summary = select_first([MosDepthOverBed.cov_stat_summary, MosDepthNoBed.cov_stat_summary])
        Array[File]? stats_per_interval_cov_stat_summary_files      = cov_stat_summary_files
        Array[Map[String, Float]]? stats_per_interval_cov_stat_summaries = cov_stat_summaries
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
        preemptible_tries: "Number of times to retry a preempted task."
    }

    input {
        File bam
        File bai
        File? bed

        String? chrom
        Int bin_length = 1000

        # Runtime parameters
        Int mem = 8
        Int preemptible_tries = 3
    }
    # mosdepth parameters
    Int threads = 4
    Boolean no_per_base = false
    Boolean fast_mode = false
    Int mapq = 1

    # coverage_stats.py parameters
    Int cov_col = 4 # column holding the coverage values
    Int round = 2
    String header_suffix = "_coverage"

    # Calculate disk size
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
        ~{prefix} ./~{basename}.bam

        # Run coverage_stats.py
        python3 /coverage_stats.py \
        --cov_col ~{cov_col} \
        --round ~{round} \
        --output_prefix ~{prefix} \
        ~{prefix}.per-base.bed.gz
    }

    output {
        # mosdepth output
        File global_dist      = "~{prefix}.mosdepth.global.dist.txt"
        File summary          = "~{prefix}.mosdepth.summary.txt"
        File? per_base        = "~{prefix}.per-base.bed.gz"
        File region_dist      = "~{prefix}.mosdepth.region.dist.txt"
        File regions          = "~{prefix}.regions.bed.gz"
        File regions_csi      = "~{prefix}.regions.bed.gz.csi"
        # coverage_stats.py output
        File cov_stat_summary_file = "~{prefix}.cov_stat_summary.json"
        Map[String, Float] cov_stat_summary = read_json("~{prefix}.cov_stat_summary.json")

    }

    runtime {
        cpu:                    4
        memory:                 mem + " GiB"
        disks: "local-disk " +  disk_size + " HDD"
        preemptible:            preemptible_tries
        maxRetries:             1
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-mosdepth:bs-cov-sum-0.3.1"
    }
}




task BinBed {
    input {
        File? bed_file
        Int? bin_length
    }

    command {
        set -euxo pipefail

        bedtools makewindows -b ~{bed_file} -w ~{bin_length} > binned.bed

    }

    output {
        File binned_bed = "binned.bed"
    }

    runtime {
        cpu: 1
        memory: "1 GiB"
        disks: "local-disk 10 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
    }
}

task FailWorkflow {
    input{
        String message
    }
    command {
        echo "Failing workflow"
        echo ~{message}
        exit 1
    }
    output {
        String out_message = message
    }
    runtime {
        docker: "marketplace.gcr.io/google/ubuntu2004"
    }
}

