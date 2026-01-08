version 1.0

import "../../structs/Structs.wdl"
import "../Utility/Finalize.wdl" as FF

workflow AlignedMetrics {

    meta {
        description: "Workflow to generate metrics for aligned BAMs"
    }
    parameter_meta {
        aligned_bam: "Aligned BAM file"
        aligned_bai: "Index for aligned BAM file"
        ref_fasta: "Reference FASTA file"
        ref_dict: "Reference dictionary file"
        gcs_output_dir: "GCS output directory"
    }

    input {
        File aligned_bam
        File aligned_bai

        File ref_fasta
        File ref_dict

        String? gcs_output_dir
    }

    call ReadMetrics as AlignedReadMetrics { input: bam = aligned_bam }

    call MakeChrIntervalList { input: ref_dict = ref_dict }

    scatter (chr_info in MakeChrIntervalList.chrs) {
        call MosDepth {
            input:
                bam = aligned_bam,
                bai = aligned_bai,
                chr = chr_info[0]
        }

        call SummarizeDepth { input: regions = MosDepth.regions }
    }

    call FlagStats as AlignedFlagStats { input: bam = aligned_bam }

    if (defined(gcs_output_dir)) {
        String outdir = sub(gcs_output_dir + "", "/$", "")

        call FF.FinalizeToDir as FFYieldAligned {
            input:
                outdir = outdir + "/yield_aligned/",
                files = [
                    AlignedFlagStats.flag_stats,
                    AlignedReadMetrics.np_hist,
                    AlignedReadMetrics.range_gap_hist,
                    AlignedReadMetrics.zmw_hist,
                    AlignedReadMetrics.prl_counts,
                    AlignedReadMetrics.prl_hist,
                    AlignedReadMetrics.prl_nx,
                    AlignedReadMetrics.prl_yield_hist,
                    AlignedReadMetrics.rl_counts,
                    AlignedReadMetrics.rl_hist,
                    AlignedReadMetrics.rl_nx,
                    AlignedReadMetrics.rl_yield_hist
                ]
        }

        call FF.FinalizeToDir as FFCoverageFullDist { input: outdir = outdir + "/coverage/", files = MosDepth.full_dist }
        call FF.FinalizeToDir as FFCoverageGlobalDist { input: outdir = outdir + "/coverage/", files = MosDepth.global_dist }
        call FF.FinalizeToDir as FFCoverageRegionDist { input: outdir = outdir + "/coverage/", files = MosDepth.region_dist }
        call FF.FinalizeToDir as FFCoverageRegions { input: outdir = outdir + "/coverage/", files = MosDepth.regions }
        call FF.FinalizeToDir as FFCoverageRegionsCsi { input: outdir = outdir + "/coverage/", files = MosDepth.regions_csi }
        call FF.FinalizeToDir as FFCoverageQuantizedDist { input: outdir = outdir + "/coverage/", files = MosDepth.quantized_dist }
        call FF.FinalizeToDir as FFCoverageQuantized { input: outdir = outdir + "/coverage/", files = MosDepth.quantized }
        call FF.FinalizeToDir as FFCoverageQuantizedCsi { input: outdir = outdir + "/coverage/", files = MosDepth.quantized_csi }

        call FF.FinalizeToDir as FFDepthSummaries { input: outdir = outdir + "/coverage_summaries/", files = SummarizeDepth.cov_summary }
    }

    output {
        File aligned_flag_stats = AlignedFlagStats.flag_stats

        Array[File] coverage_full_dist      = MosDepth.full_dist
        Array[File] coverage_global_dist    = MosDepth.global_dist
        Array[File] coverage_region_dist    = MosDepth.region_dist
        Array[File] coverage_regions        = MosDepth.regions
        Array[File] coverage_regions_csi    = MosDepth.regions_csi
        Array[File] coverage_quantized_dist = MosDepth.quantized_dist
        Array[File] coverage_quantized      = MosDepth.quantized
        Array[File] coverage_quantized_csi  = MosDepth.quantized_csi

        File aligned_np_hist = AlignedReadMetrics.np_hist
        File aligned_range_gap_hist = AlignedReadMetrics.range_gap_hist
        File aligned_zmw_hist = AlignedReadMetrics.zmw_hist
        File aligned_prl_counts = AlignedReadMetrics.prl_counts
        File aligned_prl_hist = AlignedReadMetrics.prl_hist
        File aligned_prl_nx = AlignedReadMetrics.prl_nx
        File aligned_prl_yield_hist = AlignedReadMetrics.prl_yield_hist
        File aligned_rl_counts = AlignedReadMetrics.rl_counts
        File aligned_rl_hist = AlignedReadMetrics.rl_hist
        File aligned_rl_nx = AlignedReadMetrics.rl_nx
        File aligned_rl_yield_hist = AlignedReadMetrics.rl_yield_hist

        File raw_chr_intervals = MakeChrIntervalList.raw_chrs
    }
}

task MakeChrIntervalList {
    input {
        File ref_dict

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    command <<<
        set -euxo pipefail

        grep '^@SQ' ~{ref_dict} | awk '{ print $2 "\t" 1 "\t" $3 }' | sed 's/[SL]N://g' | grep -v -e random -e chrUn -e decoy -e alt -e HLA -e EBV > chrs.txt
    >>>

    output {
        Array[Array[String]] chrs = read_tsv("chrs.txt")
        File raw_chrs = "chrs.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MosDepth {
    input {
        File bam
        File bai
        String chr
        Int? window_size

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB"))
    Int ws = select_first([window_size, 500])
    String basename = basename(bam, ".bam")
    String prefix = "~{basename}.coverage.~{chr}"

    command <<<
        set -euxo pipefail

        mosdepth -t 4 -c "~{chr}" -n -x -Q 1 ~{prefix}.full ~{bam}
        mosdepth -t 4 -c "~{chr}" -n -x -Q 1 -b ~{ws} ~{prefix} ~{bam}

        export MOSDEPTH_Q0=NO_COVERAGE   # 0 -- defined by the arguments to --quantize
        export MOSDEPTH_Q1=LOW_COVERAGE  # 1..4
        export MOSDEPTH_Q2=CALLABLE      # 5..149
        export MOSDEPTH_Q3=HIGH_COVERAGE # 150 ...

        mosdepth -t 4 -c "~{chr}" -n -x -Q 1 --quantize 0:1:5:150: ~{prefix}.quantized ~{bam}
    >>>

    output {
        File full_dist      = "~{prefix}.full.mosdepth.global.dist.txt"
        File global_dist    = "~{prefix}.mosdepth.global.dist.txt"
        File region_dist    = "~{prefix}.mosdepth.region.dist.txt"
        File regions        = "~{prefix}.regions.bed.gz"
        File regions_csi    = "~{prefix}.regions.bed.gz.csi"
        File quantized_dist = "~{prefix}.quantized.mosdepth.global.dist.txt"
        File quantized      = "~{prefix}.quantized.quantized.bed.gz"
        File quantized_csi  = "~{prefix}.quantized.quantized.bed.gz.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-mosdepth:0.3.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MosDepthOverBed {
    input {
        File bam
        File bai
        File bed

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB"))
    String basename = basename(bam, ".bam")
    String bedname = basename(bed, ".bed")
    String prefix = "~{basename}.coverage_over_bed.~{bedname}"

    command <<<
        set -euxo pipefail

        mosdepth -t 4 -b ~{bed} -n -x -Q 1 ~{prefix} ~{bam}
    >>>

    output {
        File global_dist      = "~{prefix}.mosdepth.global.dist.txt"
        File region_dist      = "~{prefix}.mosdepth.region.dist.txt"
        File regions          = "~{prefix}.regions.bed.gz"
        File regions_csi      = "~{prefix}.regions.bed.gz.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/mosdepth:0.3.3--h37c5b7d_2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MosDepthWGS {
    meta {
        description: "Collects WGS coverage of the bam. Optionally, collects coverage number for each region in the provided BED (this avoids extra localizations)."
    }
    parameter_meta {
        bed_descriptor: "A short description on the BED file provided. It will be used in naming the regions output, so be careful what you provide here."
        regions: "When bed is provided, this gets generated, which holds the coverage over the regions defined in the bed file."
    }
    input {
        File bam
        File bai
        File? bed
        String bed_descriptor = "unknown"

        String disk_type = "SSD"
        RuntimeAttr? runtime_attr_override
    }

    output {
        Float wgs_cov = read_float("wgs.cov.txt")
        File summary_txt = "~{prefix}.mosdepth.summary.txt"
        File? regions = "~{prefix}.coverage_over_bed.~{bed_descriptor}.regions.bed.gz"
    }

    String basename = basename(bam, ".bam")
    String prefix = "~{basename}.mosdepth_coverage"

    Boolean collect_over_bed = defined(bed)

    command <<<
    set -euxo pipefail

        mosdepth \
            -t 2 \
            -x -n -Q1 \
            ~{prefix} \
            ~{bam} &

        if ~{collect_over_bed}; then
            mosdepth \
                -t 2 \
                -b ~{bed} \
                -x -n -Q 1 \
                "~{prefix}.coverage_over_bed.~{bed_descriptor}" \
                ~{bam} &
        fi

        wait && ls

        # wg
        tail -n1 ~{prefix}.mosdepth.summary.txt | \
            awk -F '\t' '{print $4}' | \
            xargs printf "%0.2f\n" > wgs.cov.txt
    >>>

    #########################
    Int pd_disk_size = 10 + ceil(size(bam, "GiB"))
    Int local_disk_size = if(size(bam, "GiB")>300) then 750 else 375
    Int disk_size = if('LOCAL'==disk_type) then local_disk_size else pd_disk_size

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/mosdepth:0.3.4-gcloud"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SummarizeDepth {
    input {
        File regions

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(regions, "GB"))
    String regName = sub(basename(regions, ".regions.bed.gz"), "out.coverage.", "")

    command <<<
        set -euxo pipefail

        ((echo 'chr start stop cov_mean cov_sd cov_q1 cov_median cov_q3 cov_iqr') && \
         (zcat ~{regions} | datamash first 1 first 2 last 3 mean 4 sstdev 4 q1 4 median 4 q3 4 iqr 4)) | \
         column -t > ~{regName}.summary.txt
    >>>

    output {
        File cov_summary = "~{regName}.summary.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CoverageTrack {
    input {
        File bam
        File bai
        String chr
        String start
        String end

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB"))

    command <<<
        set -euxo pipefail

        samtools depth -a ~{bam} -r ~{chr}:~{start}-~{end} | bgzip > ~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz
        tabix -p bed ~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz
    >>>

    output {
        File coverage = "~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz"
        File coverage_tbi = "~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FlagStats {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools flagstat ~{bam} > ~{basename}.flag_stats.txt
    >>>

    output {
        File flag_stats = "~{basename}.flag_stats.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ReadNamesAndLengths {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools view ~{bam} | awk '{ print $1, length($10) }' | gzip -1 > ~{basename}.read_names_and_lengths.txt.gz
    >>>

    output {
        File read_names_and_lengths = "~{basename}.read_names_and_lengths.txt.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FilterMQ0Reads {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB"))
    String prefix = basename(bam, ".bam")

    command <<<
        set -euxo pipefail

        samtools view -q 1 -b ~{bam} > ~{prefix}.no_mq0.bam
        samtools index ~{prefix}.no_mq0.bam
    >>>

    output {
        File no_mq0_bam = "~{prefix}.no_mq0.bam"
        File no_mq0_bai = "~{prefix}.no_mq0.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ComputeBedCoverage {
    input {
        File bam
        File bai
        File bed
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB") + size(bed, "GB"))

    command <<<
        set -euxo pipefail

        bedtools coverage -b ~{bed} -a ~{bam} -nobuf | gzip > ~{prefix}.txt.gz
        zcat ~{prefix}.txt.gz | awk '{ sum += sprintf("%f", $15*$16) } END { printf("%f\n", sum) }' > ~{prefix}.count.txt
    >>>

    output {
        File coverage = "~{prefix}.txt.gz"
        Float counts = read_float("~{prefix}.count.txt")
        File counts_file = "~{prefix}.count.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ReadMetrics {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        java -jar /usr/local/bin/gatk.jar ComputeLongReadMetrics -I ~{bam} -O ~{basename}.read_metrics -DF WellformedReadFilter
    >>>

    output {
        File np_hist = "~{basename}.read_metrics.np_hist.txt"
        File range_gap_hist = "~{basename}.read_metrics.range_gap_hist.txt"
        File zmw_hist = "~{basename}.read_metrics.zmw_hist.txt"
        File prl_counts = "~{basename}.read_metrics.prl_counts.txt"
        File prl_hist = "~{basename}.read_metrics.prl_hist.txt"
        File prl_nx = "~{basename}.read_metrics.prl_nx.txt"
        File prl_yield_hist = "~{basename}.read_metrics.prl_yield_hist.txt"
        File rl_counts = "~{basename}.read_metrics.rl_counts.txt"
        File rl_hist = "~{basename}.read_metrics.rl_hist.txt"
        File rl_nx = "~{basename}.read_metrics.rl_nx.txt"
        File rl_yield_hist = "~{basename}.read_metrics.rl_yield_hist.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             50,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task BamToBed {
    input {
        File bam
        File bai

        RuntimeAttr? runtime_attr_override
    }

    String bed = basename(bam, ".bam") + ".bed"
    Int disk_size = 4*ceil(size(bam, "GB") + size(bai, "GB"))

    command <<<
        set -euxo pipefail

        bedtools bamtobed -i ~{bam} > ~{bed}
    >>>

    output {
        File bedfile = bed
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
