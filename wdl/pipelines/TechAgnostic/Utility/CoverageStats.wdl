version 1.0

workflow MosdepthCoverageStats {

    meta {
        description: "Calculate coverage statistics using mosdepth."
    }
    parameter_meta {
        aligned_bam: "Aligned BAM file."
        aligned_bai: "Aligned BAM index file."
        bed_file: "BED file containing regions of interest. Workflow will fail if bed file is empty"
        bin_length: "Length of bins to use for coverage calculation."
        preemptible_tries: "Number of times to retry a preempted task."
        stats_per_interval: "If true, calculate coverage statistics per interval. Must provide a bed file. Workflow will fail if bed file contains more than 200 intervals to avoid accidently over spending."
        max_bed_intervals: "Maximum number of intervals allowed in the bed file for MosdpethPerInterval. Default is 200."
    }

    input {
        File aligned_bam
        File aligned_bai
        File? bed_file

        Int? bin_length

        Boolean stats_per_interval = false
        Int max_bed_intervals = 200

        # Runtime parameters
        Int preemptible_tries = 3
        Int mosdepth_mem = 8
        Int mosdepth_extra_disk = 0
    }

    if (defined(bed_file)) {
        if (length(read_lines(select_first([bed_file]))) == 0) {
            call FailWorkflow as EmptyBedFile {
                input:
                    message = "stats_per_interval is set to true, but the provided bed file is empty."
            }
        }
    }

    if (stats_per_interval) {
        if (!defined(bed_file)) {
            call FailWorkflow as NoBedFile {
                input:
                    message = "stats_per_interval is set to true, but no bed file is provided."
            }
        }
        if (length(read_lines(select_first([bed_file]))) > max_bed_intervals) {
            call FailWorkflow as TooManyIntervals {
                input:
                    message = "stats_per_interval is set to true, but the provided bed file contains more than 200 intervals."
            }
        }
            call MosDepthPerInterval {
                input:
                    bam = aligned_bam,
                    bai = aligned_bai,
                    bed = bed_file,
                    preemptible_tries = preemptible_tries,
                    mem = mosdepth_mem,
                    extra_disk = mosdepth_extra_disk
            }
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
                extra_disk = mosdepth_extra_disk,
                summarize_regions = true,
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
                extra_disk = mosdepth_extra_disk,
        }
    }


    output {
        File cov_stat_summary_file      = select_first([MosDepthOverBed.cov_stat_summary_file, MosDepthNoBed.cov_stat_summary_file])
        Map[String, Float] cov_stat_summary = select_first([MosDepthOverBed.cov_stat_summary, MosDepthNoBed.cov_stat_summary])
        File? stats_per_interval_cov_stat_summary_files      = MosDepthPerInterval.cov_stat_summary_all_file
#        Array[Map[String, Object]]? stats_per_interval_cov_stat_summaries = MosDepthPerInterval.cov_stat_summary_all
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
        extra_disk: "Extra disk space to allocate for intermediate files."
    }

    input {
        File bam
        File bai
        File? bed

        String? chrom
        Int bin_length = 1000
        Boolean summarize_regions = false

        # Runtime parameters
        Int mem = 8
        Int preemptible_tries = 3
        Int extra_disk = 0
    }
    # mosdepth parameters
    Int threads = 4
    Boolean no_per_base = false
    Boolean fast_mode = false
    Int mapq = 1

    # coverage_stats.py parameters
    Int cov_col = 4 # column holding the coverage values
    Int round = 2

    # Calculate disk size
    Int disk_size = (2 * ceil(size(bam, "GB") + size(bai, "GB"))) + 100 + extra_disk
    String basename = basename(bam, ".bam")
    String prefix = "~{basename}.coverage_over_bed"

    String cov_file_to_summarize = if summarize_regions then "~{prefix}.regions.bed.gz" else "~{prefix}.per-base.bed.gz"

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
        ~{cov_file_to_summarize}
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
        memory:                 mem + " GiB"
        disks: "local-disk " +  disk_size + " HDD"
        preemptible:            preemptible_tries
        maxRetries:             1
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-mosdepth:bs-cov-sum-0.3.1"
    }
}

task MosDepthPerInterval {

    meta {
        description: "Calculate coverage using mosdepth for each interval in bed."
    }
    parameter_meta {
        bam: "Aligned BAM file."
        bai: "Aligned BAM index file."
        bed: "BED file containing regions of interest."
        bin_length: "Length of bins to use for coverage calculation. default is 1000. If bed file is provided, this parameter is ignored."
        preemptible_tries: "Number of times to retry a preempted task."
        extra_disk: "Extra disk space to allocate for intermediate files."
    }

    input {
        File bam
        File bai
        File? bed

        Int bin_length = 1000

        # Runtime parameters
        Int mem = 8
        Int preemptible_tries = 3
        Int extra_disk = 0
    }
    # mosdepth parameters
    Int threads = 4
    Boolean no_per_base = false
    Boolean fast_mode = false
    Int mapq = 1

    # coverage_stats.py parameters
    Int cov_col = 4 # column holding the coverage values
    Int round = 2

    # Calculate disk size
    Int disk_size = (2 * ceil(size(bam, "GB") + size(bai, "GB"))) + 100 + extra_disk
    String basename = basename(bam, ".bam")
    String prefix = "~{basename}.coverage_over_bed"

    String cov_file_to_summarize = "~{prefix}.per-base.bed.gz"

    command <<<
        set -euxo pipefail

        # Create symbolic links for bam and bai in the current working directory
        ln -s ~{bam} ./~{basename}.bam
        ln -s ~{bai} ./~{basename}.bai

        # Create file for coverage stats summary of all intervals
        touch ~{prefix}.cov_stat_summary_all.json
        echo "[" > ~{prefix}.cov_stat_summary_all.json

        # Create a temporary directory for intermediate files
        tmp_dir=$(mktemp -d)
        trap "rm -rf $tmp_dir" EXIT

        mapfile -t lines < "~{bed}"

        for line in "${lines[@]}"; do
            bed_file="$tmp_dir/bed_line.bed"
            echo $line | tr ' ' '\t' > $bed_file

            # Use samtools to extract reads overlapping the interval
            samtools view -bM -L $bed_file ~{basename}.bam > $tmp_dir/interval.bam
            samtools index $tmp_dir/interval.bam

            mosdepth \
            ~{true="-n" false="" no_per_base} \
            ~{true="-x" false="" fast_mode} \
            -t ~{threads} \
            ~{"-Q " + mapq} \
            ~{prefix} $tmp_dir/interval.bam

            # if mosdepth fails, exit with failure
            if [[ $? -ne 0 ]]; then
                echo "Mosdepth failed for interval $line"
                exit 1
            fi

            # Run coverage_stats.py
            python3 /coverage_stats.py \
            --debug \
            --cov_col ~{cov_col} \
            --round ~{round} \
            --output_prefix ~{prefix} \
            ~{cov_file_to_summarize}

            # In the coverage stats summary replace the open bracket with 'interval' name and the line as the value
            sed -i "s/{/\{\"interval\": \"$line\", /" ~{prefix}.cov_stat_summary.json


            # Append the coverage stats summary of the current interval to the file containing the summary of all intervals
            if [[ -s ~{prefix}.cov_stat_summary.json ]]; then
                cat ~{prefix}.cov_stat_summary.json >> ~{prefix}.cov_stat_summary_all.json
                echo "," >> ~{prefix}.cov_stat_summary_all.json
            fi

        done
        # remove the last comma and close the json array
        sed -i '$ s/,$//' ~{prefix}.cov_stat_summary_all.json
        echo "]" >> ~{prefix}.cov_stat_summary_all.json


    >>>

    output {
        # coverage_stats.py output
        File cov_stat_summary_all_file = "~{prefix}.cov_stat_summary_all.json"
#        Array[Map[String, Object]] cov_stat_summary_all = read_json("~{prefix}.cov_stat_summary_all.json")

    }

    runtime {
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
        set -e

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

