version 1.0

import "../../structs/Structs.wdl"

workflow AlignedMetrics {

    meta {
        description: "Workflow to generate metrics for aligned BAMs"
    }
    parameter_meta {
        aligned_bam: "Aligned BAM file"
        aligned_bai: "Index for aligned BAM file"
        ref_fasta: "Not used"
        ref_dict: "Reference dictionary file"
        gcs_output_dir: "GCS output directory"
    }

    input {
        File aligned_bam
        File aligned_bai

        File? ref_fasta
        File ref_dict

        String? gcs_output_dir
    }

    call ReadMetrics {
        input:
            bam = aligned_bam,
            gcs_output_dir = gcs_output_dir
    }

    call MakeChrIntervalList { input: ref_dict = ref_dict }

    scatter (chr_info in MakeChrIntervalList.chrs) {
        call MosDepth {
            input:
                bam = aligned_bam,
                bai = aligned_bai,
                chr = chr_info[0],
                gcs_output_dir = gcs_output_dir
        }
    }

    output {
        Array[File] coverage_full_dist      = MosDepth.full_dist
        Array[File] coverage_global_dist    = MosDepth.global_dist
        Array[File] coverage_region_dist    = MosDepth.region_dist
        Array[File] coverage_regions        = MosDepth.regions
        Array[File] coverage_regions_csi    = MosDepth.regions_csi
        Array[File] coverage_quantized_dist = MosDepth.quantized_dist
        Array[File] coverage_quantized      = MosDepth.quantized
        Array[File] coverage_quantized_csi  = MosDepth.quantized_csi
        Array[File] coverage_summary        = MosDepth.cov_summary

        File aligned_np_hist = ReadMetrics.np_hist
        File aligned_range_gap_hist = ReadMetrics.range_gap_hist
        File aligned_zmw_hist = ReadMetrics.zmw_hist
        File aligned_prl_counts = ReadMetrics.prl_counts
        File aligned_prl_hist = ReadMetrics.prl_hist
        File aligned_prl_nx = ReadMetrics.prl_nx
        File aligned_prl_yield_hist = ReadMetrics.prl_yield_hist
        File aligned_rl_counts = ReadMetrics.rl_counts
        File aligned_rl_hist = ReadMetrics.rl_hist
        File aligned_rl_nx = ReadMetrics.rl_nx
        File aligned_rl_yield_hist = ReadMetrics.rl_yield_hist
        File aligned_flag_stats = ReadMetrics.flag_stats

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
        docker:             "ubuntu:21.04"
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
        Int window_size = 500
        String? gcs_output_dir

        RuntimeAttr? runtime_attr_override
    }

    String? coverage_dir = if defined(gcs_output_dir) then sub(select_first([gcs_output_dir]), "/?$", "/coverage/") else gcs_output_dir
    String? summary_dir = if defined(gcs_output_dir) then sub(select_first([gcs_output_dir]), "/?$", "/coverage_summaries/") else gcs_output_dir
    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB"))
    String basename = basename(bam, ".bam")
    String prefix = "~{basename}.coverage.~{chr}"

    command <<<
        set -euxo pipefail

        mosdepth -t 4 -c "~{chr}" -n -x -Q 1 ~{prefix}.full ~{bam}
        mosdepth -t 4 -c "~{chr}" -n -x -Q 1 -b ~{window_size} ~{prefix} ~{bam}

        export MOSDEPTH_Q0=NO_COVERAGE   # 0 -- defined by the arguments to --quantize
        export MOSDEPTH_Q1=LOW_COVERAGE  # 1..4
        export MOSDEPTH_Q2=CALLABLE      # 5..149
        export MOSDEPTH_Q3=HIGH_COVERAGE # 150 ...

        mosdepth -t 4 -c "~{chr}" -n -x -Q 1 --quantize 0:1:5:150: ~{prefix}.quantized ~{bam}

        ( echo 'chr start stop cov_mean cov_sd cov_q1 cov_median cov_q3 cov_iqr' && \
          zcat ~{prefix}.regions.bed.gz | datamash first 1 first 2 last 3 mean 4 sstdev 4 q1 4 median 4 q3 4 iqr 4 ) | \
            column -t > ~{prefix}.summary.txt

        if ~{defined(coverage_dir)}; then
          gsutil -m cp \
             "~{prefix}.full.mosdepth.global.dist.txt" \
             "~{prefix}.mosdepth.global.dist.txt" \
             "~{prefix}.mosdepth.region.dist.txt" \
             "~{prefix}.regions.bed.gz" \
             "~{prefix}.regions.bed.gz.csi" \
             "~{prefix}.quantized.mosdepth.global.dist.txt" \
             "~{prefix}.quantized.quantized.bed.gz" \
             "~{prefix}.quantized.quantized.bed.gz.csi" \
             "~{coverage_dir}"
          gsutil cp "~{prefix}.summary.txt" "~{summary_dir}"
        fi

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
        File cov_summary    = "~{prefix}.summary.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-mosdepth:0.3.2"
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
        String? gcs_output_dir

        RuntimeAttr? runtime_attr_override
    }

    String? output_dir = if defined(gcs_output_dir) then sub(select_first([gcs_output_dir]), "/?$", "/yield_aligned/") else gcs_output_dir
    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        cat << "EOF" > summarizeBAM.awk
        BEGIN { FS=OFS="\t"; lastZMWid =-1; }
        {
                # retrieve the number of passes as the value of the "np" tag (default = 1)
                nPasses=1;
                for( fld=12;fld<=NF;++fld )
                        if( substr($fld,1,5)=="np:i:" ) {
                                nPasses=substr($fld,6);
                                break;
                        }

                # get the read length and its bin for histograms
                readLength=length($10);

                # default the zero-mode waveguide id as last value + 1
                zmwId=lastZMWid+1;

                #default the range info as 0:0
                rangeStart=0;
                rangeStop=0;

                # parse the read name
                if( split($1,namePieces,"/")==3 ) {
                        zmwId = namePieces[2];
                        if( split(namePieces[3],range,"_")==2 ) {
                                rangeStart = range[0];
                                rangeStop = range[1];
                        }
                }

                # calculate value to report for range gap (or -1 to indicate no report for this read)
                reportedRangeGap=-1;
                if( lastRangeStop!=0 && lastZMWid==zmwId ) {
                        reportedRangeGap=rangeStart-lastRangeStop;
                }

                # calculate value to report for polymerase read length (or -1 to indicate no report for this read)
                reportedPRL = -1;
                if( lastZMWid!=-1 && lastZMWid!=zmwId ) {
                        reportedPRL=polymeraseRL;
                }

                # track polymerase read length across reads with the same zmwId
                if( zmwId!=lastZMWid ) polymeraseRL=0;
                polymeraseRL+=readLength;

                # update values maintained across reads
                lastZMWid=zmwId;
                lastRangeStop=rangeStop;

                # print the parsed data
                print readLength,nPasses,readLength*nPasses,zmwId,reportedRangeGap,reportedPRL;
        }
        EOF

        cat << "EOF" > calculateNX.awk
        BEGIN { FS=OFS="\t" }
        {
                sum+=$1;
                cent=int(100.*sum/tot);
                if( cent!=lastCent )
                        print cent,lastVal;
                lastCent=cent;
                lastVal=$1
        }
        EOF

        # awkOut.tsv file is a tab delimited file with these columns:
        # fld:       1        2             3           4          5                      6
        # datum: readLength nPasses yield(len*passes) zmwId reportedRangeGap reportedPolymeraseReadLength
        samtools view -F0x900 ~{bam} | awk -f summarizeBAM.awk > awkOut.tsv

        # extract the nPasses column for each read, sort numerically, count each group, add a header line
        cut -f2 awkOut.tsv | sort -n | datamash -g1 count 1 | \
            awk 'BEGIN{print "np\tcount"}{print}' > "~{basename}.read_metrics.np_hist.txt"

        # extract the reportedRangeGap column, ignore -1 values, sort numerically, count each group, add header line
        cut -f5 awkOut.tsv | sed '/^-1$/d' | sort -n | datamash -g1 count 1 | \
            awk 'BEGIN{print "range_gap\tcount"}{print}' > "~{basename}.read_metrics.range_gap_hist.txt"

        # extract the zmwId column, sort numerically, count each group, extract counts and sort numerically,
        #    count each group of counts, add header line
        cut -f4 awkOut.tsv | sort -n | datamash -g1 count 1 | cut -f2 | sort -n | datamash -g1 count 1 | \
            awk 'BEGIN{print "zmw_count\tcount"}{print}' > "~{basename}.read_metrics.zmw_hist.txt"

        # extract yield and prl, ignore -1 prls, sort numerically, mash some summary stats, add header line
        cut -f3,6 awkOut.tsv | sed '/	-1$/d' | sort -n | datamash count 2 sum 2 sum 1 min 2 max 2 mean 2 sstdev 2 | \
            awk 'BEGIN{print "num_polymerase_reads\tpolymerase_read_yield_bp\tpolymerase_read_yield_bp_times_passes\tpolymerase_read_min_length_bp\tpolymerase_read_max_length_bp\tpolymerase_read_mean_length_bp\tpolymerase_read_sd_length_bp"}{print}' > "~{basename}.read_metrics.prl_counts.txt"

        # extract prl rounded to 100s ignoring -1 values, sort numerically, count each group, add header line
        awk 'BEGIN{FS=OFS="\t"}{if($6>=0)print int(($6+99)/100)*100}' awkOut.tsv | sort -n | datamash -g1 count 1 | \
            awk 'BEGIN{print "polymerase_read_length\tlength_bp"}{print}' > "~{basename}.read_metrics.prl_hist.txt"

        # calculate total of prl column, ignoring -1 values
        # then extract the prl column, ignore -1 values, reverse sort numerically, print the centile divisions,
        #    forget about the 100th percentile, add header line
        tot=$(cut -f6 < awkOut.tsv | sed '/^-1$/d' | datamash sum 1)
        cut -f6 < awkOut.tsv | sed '/^-1$/d' | sort -rn | awk -v "tot=$tot" -f calculateNX.awk | sed '/^100	/d' | \
            awk 'BEGIN{print "NX\tvalue"}{print}' > "~{basename}.read_metrics.prl_nx.txt"

        # extract prl rounded to 100s and the original prl, sort numerically, sum prls for each bin, add header line
        awk 'BEGIN{FS=OFS="\t"}{if($6>=0)print int(($6+99)/100)*100,$6}' awkOut.tsv | sort -n | datamash -g1 sum 2 | \
            awk 'BEGIN{print "polymerase_read_length\tyield_bp"}{print}' > "~{basename}.read_metrics.prl_yield_hist.txt"

        # extract readLength and yield columns, sort numerically, mash some summary stats, add header line
        cut -f1,3 awkOut.tsv | sort -n | datamash count 1 sum 1 sum 2 min 1 max 1 mean 1 sstdev 1 | \
            awk 'BEGIN{print "num_reads\tread_yield_bp\tread_yield_bp_times_passes\tread_min_length_bp\tread_max_length_bp\tread_mean_length_bp\tread_sd_length_bp"}{print}' > "~{basename}.read_metrics.rl_counts.txt"

        # extract readLength rounded to 100s, sort numerically, count each group, add header line
        awk 'BEGIN{FS=OFS="\t"}{print int(($1+99)/100)*100}' awkOut.tsv | sort -n | datamash -g1 count 1 | \
            awk 'BEGIN{print "read_length_bp\tcount"}{print}' > "~{basename}.read_metrics.rl_hist.txt"

        # calculate total of readLength column,
        # then extract the readLength column, reverse sort numerically, print the centile divisions,
        #    forget about the 100th percentile, add header line
        tot=$(cut -f1 < awkOut.tsv | datamash sum 1)
        cut -f1 < awkOut.tsv | sort -rn | awk -v "tot=$tot" -f calculateNX.awk | sed '/^100	/d' | \
            awk 'BEGIN{print "NX\tvalue"}{print}' > "~{basename}.read_metrics.rl_nx.txt"

        # extract readLength rounded to 100s and the original readLength, sort numerically, sum each group, add header line
        awk 'BEGIN{FS=OFS="\t"}{print int(($1+99)/100)*100,$1}' awkOut.tsv | sort -n | datamash -g1 sum 2 | \
            awk 'BEGIN{print "read_length\tyield_bp"}{print}' > "~{basename}.read_metrics.rl_yield_hist.txt"

        samtools flagstat ~{bam} > ~{basename}.flag_stats.txt

        if ~{defined(output_dir)}; then
          gsutil -m cp \
            "~{basename}.read_metrics.np_hist.txt" \
            "~{basename}.read_metrics.range_gap_hist.txt" \
            "~{basename}.read_metrics.zmw_hist.txt" \
            "~{basename}.read_metrics.prl_counts.txt" \
            "~{basename}.read_metrics.prl_hist.txt" \
            "~{basename}.read_metrics.prl_nx.txt" \
            "~{basename}.read_metrics.prl_yield_hist.txt" \
            "~{basename}.read_metrics.rl_counts.txt" \
            "~{basename}.read_metrics.rl_hist.txt" \
            "~{basename}.read_metrics.rl_nx.txt" \
            "~{basename}.read_metrics.rl_yield_hist.txt" \
            "~{basename}.flag_stats.txt" \
            "~{output_dir}"
        fi
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
        File flag_stats = "~{basename}.flag_stats.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             50,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-mini:0.1.1"
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
