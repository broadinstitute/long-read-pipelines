version 1.0

import "Structs.wdl"
import "Finalize.wdl" as FF

workflow AlignedMetrics {
    input {
        File aligned_bam
        File aligned_bai

        File ref_fasta
        File ref_dict
        File ref_flat

        File dbsnp_vcf
        File dbsnp_tbi
        File metrics_locus

        String per
        String type
        String label

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

    call ReadNamesAndLengths { input: bam = aligned_bam }

    call RnaSeqMetrics {
        input:
            bam = aligned_bam,
            bai = aligned_bai,
            ref_flat = ref_flat
    }

    call FilterMQ0Reads { input: bam = aligned_bam }

#    call CollectSamErrorMetrics {
#        input:
#            bam = aligned_bam,
#            bai = aligned_bai,
#            ref_fasta = ref_fasta,
#            ref_dict = ref_dict,
#            dbsnp_vcf = dbsnp_vcf,
#            dbsnp_tbi = dbsnp_tbi,
#            region = metrics_locus,
#    }

    if (defined(gcs_output_dir)) {
        String outdir = sub(gcs_output_dir + "", "/$", "") + "/metrics/~{per}_~{type}_~{label}"

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

        call FF.FinalizeToDir as FFRnaSeqMetrics { input: outdir = outdir + "/rnaseq/", files = [ RnaSeqMetrics.rna_metrics ] }
        call FF.FinalizeToDir as FFReadNamesAndLengths { input: outdir = outdir + "/read_names_and_lengths/", files = [ ReadNamesAndLengths.read_names_and_lengths ] }

#        call FF.FinalizeToDir as FFErrorStats {
#            input:
#                outdir = outdir + "/error_rate/",
#                files = [
#                    CollectSamErrorMetrics.error_by_all,
#                    CollectSamErrorMetrics.error_by_base_quality,
#                    CollectSamErrorMetrics.error_by_binned_length_homopolymer_and_following_ref_base,
#                    CollectSamErrorMetrics.error_by_cycle,
#                    CollectSamErrorMetrics.error_by_gc,
#                    CollectSamErrorMetrics.error_by_homopolymer_and_following_ref_base,
#                    CollectSamErrorMetrics.error_by_insert_length,
#                    CollectSamErrorMetrics.error_by_mapping_quality,
#                    CollectSamErrorMetrics.error_by_mismatches_in_read,
#                    CollectSamErrorMetrics.error_by_one_base_padded_context,
#                    CollectSamErrorMetrics.error_by_pair_orientation,
#                    CollectSamErrorMetrics.error_by_read_direction,
#                    CollectSamErrorMetrics.error_by_read_group,
#                    CollectSamErrorMetrics.error_by_read_ordinality,
#                    CollectSamErrorMetrics.error_by_read_ordinality_and_cycle,
#                    CollectSamErrorMetrics.error_by_read_ordinality_and_gc,
#                    CollectSamErrorMetrics.error_by_read_ordinality_and_homopolymer_and_following_ref_base,
#                    CollectSamErrorMetrics.error_by_read_ordinality_and_pre_dinuc,
#                    CollectSamErrorMetrics.indel_error_by_all,
#                    CollectSamErrorMetrics.overlapping_error_by_all,
#                    CollectSamErrorMetrics.overlapping_error_by_base_quality,
#                    CollectSamErrorMetrics.overlapping_error_by_insert_length,
#                    CollectSamErrorMetrics.overlapping_error_by_read_ordinality,
#                    CollectSamErrorMetrics.overlapping_error_by_read_ordinality_and_cycle,
#                    CollectSamErrorMetrics.overlapping_error_by_read_ordinality_and_gc,
#                    CollectSamErrorMetrics.overlapping_error_by_read_ordinality_and_homopolymer_and_following_ref_base
#                ]
#        }
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

        File rna_metrics = RnaSeqMetrics.rna_metrics

#        File error_by_all = CollectSamErrorMetrics.error_by_all
#        File error_by_base_quality = CollectSamErrorMetrics.error_by_base_quality
#        File error_by_binned_length_homopolymer_and_following_ref_base = CollectSamErrorMetrics.error_by_binned_length_homopolymer_and_following_ref_base
#        File error_by_cycle = CollectSamErrorMetrics.error_by_cycle
#        File error_by_gc = CollectSamErrorMetrics.error_by_gc
#        File error_by_homopolymer_and_following_ref_base = CollectSamErrorMetrics.error_by_homopolymer_and_following_ref_base
#        File error_by_insert_length = CollectSamErrorMetrics.error_by_insert_length
#        File error_by_mapping_quality = CollectSamErrorMetrics.error_by_mapping_quality
#        File error_by_mismatches_in_read = CollectSamErrorMetrics.error_by_mismatches_in_read
#        File error_by_one_base_padded_context = CollectSamErrorMetrics.error_by_one_base_padded_context
#        File error_by_pair_orientation = CollectSamErrorMetrics.error_by_pair_orientation
#        File error_by_read_direction = CollectSamErrorMetrics.error_by_read_direction
#        File error_by_read_group = CollectSamErrorMetrics.error_by_read_group
#        File error_by_read_ordinality = CollectSamErrorMetrics.error_by_read_ordinality
#        File error_by_read_ordinality_and_cycle = CollectSamErrorMetrics.error_by_read_ordinality_and_cycle
#        File error_by_read_ordinality_and_gc = CollectSamErrorMetrics.error_by_read_ordinality_and_gc
#        File error_by_read_ordinality_and_homopolymer_and_following_ref_base = CollectSamErrorMetrics.error_by_read_ordinality_and_homopolymer_and_following_ref_base
#        File error_by_read_ordinality_and_pre_dinuc = CollectSamErrorMetrics.error_by_read_ordinality_and_pre_dinuc
#        File indel_error_by_all = CollectSamErrorMetrics.indel_error_by_all
#        File overlapping_error_by_all = CollectSamErrorMetrics.overlapping_error_by_all
#        File overlapping_error_by_base_quality = CollectSamErrorMetrics.overlapping_error_by_base_quality
#        File overlapping_error_by_insert_length = CollectSamErrorMetrics.overlapping_error_by_insert_length
#        File overlapping_error_by_read_ordinality = CollectSamErrorMetrics.overlapping_error_by_read_ordinality
#        File overlapping_error_by_read_ordinality_and_cycle = CollectSamErrorMetrics.overlapping_error_by_read_ordinality_and_cycle
#        File overlapping_error_by_read_ordinality_and_gc = CollectSamErrorMetrics.overlapping_error_by_read_ordinality_and_gc
#        File overlapping_error_by_read_ordinality_and_homopolymer_and_following_ref_base = CollectSamErrorMetrics.overlapping_error_by_read_ordinality_and_homopolymer_and_following_ref_base
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
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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

        mosdepth -t 4 -c ~{chr} -n -x -Q 1 ~{prefix}.full ~{bam}
        mosdepth -t 4 -c ~{chr} -n -x -Q 1 -b ~{ws} ~{prefix} ~{bam}

        export MOSDEPTH_Q0=NO_COVERAGE   # 0 -- defined by the arguments to --quantize
        export MOSDEPTH_Q1=LOW_COVERAGE  # 1..4
        export MOSDEPTH_Q2=CALLABLE      # 5..149
        export MOSDEPTH_Q3=HIGH_COVERAGE # 150 ...

        mosdepth -t 4 -c ~{chr} -n -x -Q 1 --quantize 0:1:5:150: ~{prefix}.quantized ~{bam}

        ls -lah
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
        docker:             "quay.io/biocontainers/mosdepth:0.2.4--he527e40_0"
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
        docker:             "quay.io/biocontainers/mosdepth:0.2.4--he527e40_0"
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

task SummarizeDepth {
    input {
        File regions

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(regions, "GB"))
    String chrName = sub(basename(regions, ".regions.bed.gz"), "out.coverage.", "")

    command <<<
        set -euxo pipefail

        ((echo 'chr start stop cov_mean cov_sd cov_q1 cov_median cov_q3 cov_iqr') && \
         (zcat ~{regions} | datamash first 1 first 2 last 3 mean 4 sstdev 4 q1 4 median 4 q3 4 iqr 4)) | \
         column -t > ~{chrName}.summary.txt
    >>>

    output {
        File cov_summary = "~{chrName}.summary.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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

task RnaSeqMetrics {
    input {
        File bam
        File bai
        File ref_flat

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB") + size(ref_flat, "GB"))

    command <<<
        set -euxo pipefail

        java -jar /usr/local/bin/picard.jar CollectRnaSeqMetrics \
            I=~{bam} \
            REF_FLAT=~{ref_flat} \
            STRAND=NONE \
            VALIDATION_STRINGENCY=LENIENT \
            O=~{basename}.rna_metrics.txt

        sed -i 1,5d ~{basename}.rna_metrics.txt
    >>>

    output {
        File rna_metrics = "~{basename}.rna_metrics.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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

task CollectSamErrorMetrics {
    input {
        File bam
        File bai
        File ref_fasta
        File ref_dict
        File dbsnp_vcf
        File dbsnp_tbi
        File region

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bam, "GB") + size(bai, "GB") + size(ref_fasta, "GB") + size(dbsnp_vcf, "GB"))

    command <<<
        set -euxo pipefail

        java -jar /usr/local/bin/picard.jar CollectSamErrorMetrics R=~{ref_fasta} I=~{bam} V=~{dbsnp_vcf} O=csem L=~{region}

        find . -name 'csem.*' -exec sed -i 1,5d '{}' \;
    >>>

    output {
        File error_by_all = "csem.error_by_all"
        File error_by_base_quality = "csem.error_by_base_quality"
        File error_by_binned_length_homopolymer_and_following_ref_base = "csem.error_by_binned_length_homopolymer_and_following_ref_base"
        File error_by_cycle = "csem.error_by_cycle"
        File error_by_gc = "csem.error_by_gc"
        File error_by_homopolymer_and_following_ref_base = "csem.error_by_homopolymer_and_following_ref_base"
        File error_by_insert_length = "csem.error_by_insert_length"
        File error_by_mapping_quality = "csem.error_by_mapping_quality"
        File error_by_mismatches_in_read = "csem.error_by_mismatches_in_read"
        File error_by_one_base_padded_context = "csem.error_by_one_base_padded_context"
        File error_by_pair_orientation = "csem.error_by_pair_orientation"
        File error_by_read_direction = "csem.error_by_read_direction"
        File error_by_read_group = "csem.error_by_read_group"
        File error_by_read_ordinality = "csem.error_by_read_ordinality"
        File error_by_read_ordinality_and_cycle = "csem.error_by_read_ordinality_and_cycle"
        File error_by_read_ordinality_and_gc = "csem.error_by_read_ordinality_and_gc"
        File error_by_read_ordinality_and_homopolymer_and_following_ref_base = "csem.error_by_read_ordinality_and_homopolymer_and_following_ref_base"
        File error_by_read_ordinality_and_pre_dinuc = "csem.error_by_read_ordinality_and_pre_dinuc"
        File indel_error_by_all = "csem.indel_error_by_all"
        File overlapping_error_by_all = "csem.overlapping_error_by_all"
        File overlapping_error_by_base_quality = "csem.overlapping_error_by_base_quality"
        File overlapping_error_by_insert_length = "csem.overlapping_error_by_insert_length"
        File overlapping_error_by_read_ordinality = "csem.overlapping_error_by_read_ordinality"
        File overlapping_error_by_read_ordinality_and_cycle = "csem.overlapping_error_by_read_ordinality_and_cycle"
        File overlapping_error_by_read_ordinality_and_gc = "csem.overlapping_error_by_read_ordinality_and_gc"
        File overlapping_error_by_read_ordinality_and_homopolymer_and_following_ref_base = "csem.overlapping_error_by_read_ordinality_and_homopolymer_and_following_ref_base"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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

task SamtoolsStats {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(6 * size(bam, "GiB"))

    String raw_stats_file = "raw_stats.txt"
    String summary_stats_file = "summary_stats.txt"
    String first_frag_qual_file = "first_fragment_quality_stats.txt"
    String last_frag_qual_file = "last_fragment_quality_stats.txt"
    String first_frag_gc_content_file = "first_fragment_gc_content_stats.txt"
    String last_frag_gc_content_file = "last_fragment_gc_content_stats.txt"
    String acgt_content_per_cycle_file = "acgt_content_per_cycle_stats.txt"
    String insert_size_file = "insert_size_stats.txt"
    String read_length_dist_file = "read_length_distribution.txt"
    String indel_distribution_file = "indel_distribution.txt"
    String indels_per_cycle_file = "indels_per_cycle_stats.txt"
    String coverage_distribution_file = "coverage_distribution.txt"
    String gc_depth_file = "gc_depth_stats.txt"

    command <<<
        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        samtools stats -@$np ~{bam} > ~{raw_stats_file}
        grep ^SN ~{raw_stats_file} | cut -f 2-  > ~{summary_stats_file}
        grep ^FFQ ~{raw_stats_file} | cut -f 2- > ~{first_frag_qual_file}
        grep ^LFQ ~{raw_stats_file} | cut -f 2- > ~{last_frag_qual_file}
        grep ^GCF ~{raw_stats_file} | cut -f 2- > ~{first_frag_gc_content_file}
        grep ^GCL ~{raw_stats_file} | cut -f 2- > ~{last_frag_gc_content_file}
        grep ^GCC ~{raw_stats_file} | cut -f 2- > ~{acgt_content_per_cycle_file}
        grep ^IS ~{raw_stats_file} | cut -f 2-  > ~{insert_size_file}
        grep ^RL ~{raw_stats_file} | cut -f 2-  > ~{read_length_dist_file}
        grep ^ID ~{raw_stats_file} | cut -f 2-  > ~{indel_distribution_file}
        grep ^IC ~{raw_stats_file} | cut -f 2-  > ~{indels_per_cycle_file}
        grep ^COV ~{raw_stats_file} | cut -f 2- > ~{coverage_distribution_file}
        grep ^GCD ~{raw_stats_file} | cut -f 2- > ~{gc_depth_file}

    >>>

    output {
        File raw_stats = raw_stats_file
        File summary_stats = summary_stats_file
        File first_frag_qual = first_frag_qual_file
        File last_frag_qual = last_frag_qual_file
        File first_frag_gc_content = first_frag_gc_content_file
        File last_frag_gc_content = last_frag_gc_content_file
        File acgt_content_per_cycle = acgt_content_per_cycle_file
        File insert_size = insert_size_file
        File read_length_dist = read_length_dist_file
        File indel_distribution = indel_distribution_file
        File indels_per_cycle = indels_per_cycle_file
        File coverage_distribution = coverage_distribution_file
        File gc_depth = gc_depth_file
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
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
