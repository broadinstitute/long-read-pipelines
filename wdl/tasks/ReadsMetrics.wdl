version 1.0

import "AlignedMetrics.wdl" as AM
import "Finalize.wdl" as FF

workflow CalculateAndFinalizeReadMetrics {
    input {
        File bam_file
        File bam_index
        File ref_dict

        String base_metrics_out_dir
    }

    String read_metrics_dir = sub(base_metrics_out_dir, "/$", "") + "/read_metrics"
    String depth_metrics_dir = sub(base_metrics_out_dir, "/$", "") + "/depth_metrics"

    # Calculate metrics that require reads and reference alone:

    # Read Metrics
    call AM.ReadMetrics { input: bam = bam_file }
    call FF.FinalizeToDir as FinalizeReadMetrics {
        input:
            outdir = read_metrics_dir,
            files = [ ReadMetrics.np_hist, ReadMetrics.range_gap_hist, ReadMetrics.zmw_hist, ReadMetrics.prl_counts, ReadMetrics.prl_hist, ReadMetrics.prl_nx, ReadMetrics.prl_yield_hist, ReadMetrics.rl_counts, ReadMetrics.rl_hist, ReadMetrics.rl_nx, ReadMetrics.rl_yield_hist ]
    }

    # Depth Stats
    call AM.MakeChrIntervalList { input: ref_dict = ref_dict }
    scatter (chr_info in MakeChrIntervalList.chrs) {
        call AM.MosDepth {
            input:
                bam = bam_file,
                bai = bam_index,
                chr = chr_info[0]
        }
        call AM.SummarizeDepth { input: regions = MosDepth.regions }

        call FF.FinalizeToDir as FinalizeDepthMetrics {
            input:
                outdir = depth_metrics_dir,
                files = [ MosDepth.full_dist, MosDepth.global_dist, MosDepth.region_dist, MosDepth.regions, MosDepth.regions_csi, MosDepth.quantized_dist, MosDepth.quantized, MosDepth.quantized_csi, SummarizeDepth.cov_summary ]
        }
    }

    # Flag Stats
    call AM.FlagStats { input: bam = bam_file }
    call FF.FinalizeToDir as FinalizeFlagStats {
        input:
            outdir = base_metrics_out_dir,
            files = [ FlagStats.flag_stats ]
    }

    # Sam Stats
    call AM.SamtoolsStats { input: bam = bam_file }
    call FF.FinalizeToDir as FinalizeSamStats {
        input:
            outdir = read_metrics_dir,
            files = [
                SamtoolsStats.raw_stats,
                SamtoolsStats.summary_stats,
                SamtoolsStats.first_frag_qual,
                SamtoolsStats.last_frag_qual,
                SamtoolsStats.first_frag_gc_content,
                SamtoolsStats.last_frag_gc_content,
                SamtoolsStats.acgt_content_per_cycle,
                SamtoolsStats.insert_size,
                SamtoolsStats.read_length_dist,
                SamtoolsStats.indel_distribution,
                SamtoolsStats.indels_per_cycle,
                SamtoolsStats.coverage_distribution,
                SamtoolsStats.gc_depth
            ]
    }

    # Read Names And Lengths
    call AM.ReadNamesAndLengths { input: bam = bam_file }
    call FF.FinalizeToDir as FinalizeReadNamesAndLengths {
        input:
            outdir = base_metrics_out_dir,
            files = [ ReadNamesAndLengths.read_names_and_lengths ]
    }

    output {
        String read_metrics_path = read_metrics_dir
        String depth_metrics_path = depth_metrics_dir

        File read_metrics_np_hist = ReadMetrics.np_hist
        File read_metrics_range_gap_hist = ReadMetrics.range_gap_hist
        File read_metrics_zmw_hist = ReadMetrics.zmw_hist
        File read_metrics_prl_counts = ReadMetrics.prl_counts
        File read_metrics_prl_hist = ReadMetrics.prl_hist
        File read_metrics_prl_nx = ReadMetrics.prl_nx
        File read_metrics_prl_yield_hist = ReadMetrics.prl_yield_hist
        File read_metrics_rl_counts = ReadMetrics.rl_counts
        File read_metrics_rl_hist = ReadMetrics.rl_hist
        File read_metrics_rl_nx = ReadMetrics.rl_nx
        File read_metrics_rl_yield_hist = ReadMetrics.rl_yield_hist

        Array[Array[String]] chrs = MakeChrIntervalList.chrs

        Array[File] full_dist      = MosDepth.full_dist
        Array[File] global_dist    = MosDepth.global_dist
        Array[File] region_dist    = MosDepth.region_dist
        Array[File] regions        = MosDepth.regions
        Array[File] regions_csi    = MosDepth.regions_csi
        Array[File] quantized_dist = MosDepth.quantized_dist
        Array[File] quantized      = MosDepth.quantized
        Array[File] quantized_csi  = MosDepth.quantized_csi

        Array[File] coverage_summary = SummarizeDepth.cov_summary

        File flag_stats = FlagStats.flag_stats

        File sam_stats_raw_stats = SamtoolsStats.raw_stats
        File sam_stats_summary_stats = SamtoolsStats.summary_stats
        File sam_stats_first_frag_qual_stats = SamtoolsStats.first_frag_qual
        File sam_stats_last_frag_qual_stats = SamtoolsStats.last_frag_qual
        File sam_stats_first_frag_gc_content_stats = SamtoolsStats.first_frag_gc_content
        File sam_stats_last_frag_gc_content_stats = SamtoolsStats.last_frag_gc_content
        File sam_stats_agct_content_per_cycle_stats = SamtoolsStats.acgt_content_per_cycle
        File sam_stats_insert_size_stats = SamtoolsStats.insert_size
        File sam_stats_read_length_dist_stats = SamtoolsStats.read_length_dist
        File sam_stats_indel_dist_stats = SamtoolsStats.indel_distribution
        File sam_stats_indels_per_cycle_stats = SamtoolsStats.indels_per_cycle
        File sam_stats_coverage_dist_stats = SamtoolsStats.coverage_distribution
        File sam_stats_gc_depth_stats = SamtoolsStats.gc_depth

        File read_names_and_lengths = ReadNamesAndLengths.read_names_and_lengths
    }
}

workflow CalculateAndFinalizeAlternateReadMetrics {
    input {
        File bam_file
        File bam_index
        File ref_dict

        String base_metrics_out_dir
    }

    String read_metrics_dir = sub(base_metrics_out_dir, "/$", "") + "/read_metrics"
    String depth_metrics_dir = sub(base_metrics_out_dir, "/$", "") + "/depth_metrics"

    # Calculate metrics that require reads and reference alone:

    # Sam Stats
    call AM.SamtoolsStats { input: bam = bam_file }
    call FF.FinalizeToDir as FinalizeSamStats {
        input:
            outdir = read_metrics_dir,
            files = [
                SamtoolsStats.raw_stats,
                SamtoolsStats.summary_stats,
                SamtoolsStats.first_frag_qual,
                SamtoolsStats.last_frag_qual,
                SamtoolsStats.first_frag_gc_content,
                SamtoolsStats.last_frag_gc_content,
                SamtoolsStats.acgt_content_per_cycle,
                SamtoolsStats.insert_size,
                SamtoolsStats.read_length_dist,
                SamtoolsStats.indel_distribution,
                SamtoolsStats.indels_per_cycle,
                SamtoolsStats.coverage_distribution,
                SamtoolsStats.gc_depth
            ]
    }

    # Commented out for now because mos-depth doesn't handle contigs with dashes in their names
    # https://github.com/brentp/mosdepth/issues/115
#    # Depth Stats
#    call AM.MakeChrIntervalList { input: ref_dict = ref_dict }
#    scatter (chr_info in MakeChrIntervalList.chrs) {
#        call AM.MosDepth {
#            input:
#                bam = bam_file,
#                bai = bam_index,
#                chr = chr_info[0]
#        }
#        call AM.SummarizeDepth { input: regions = MosDepth.regions }
#
#        call FF.FinalizeToDir as FinalizeDepthMetrics {
#            input:
#                outdir = depth_metrics_dir,
#                files = [ MosDepth.full_dist, MosDepth.global_dist, MosDepth.region_dist, MosDepth.regions, MosDepth.regions_csi, MosDepth.quantized_dist, MosDepth.quantized, MosDepth.quantized_csi, SummarizeDepth.cov_summary ]
#        }
#    }

    # Flag Stats
    call AM.FlagStats { input: bam = bam_file }
    call FF.FinalizeToDir as FinalizeFlagStats {
        input:
            outdir = base_metrics_out_dir,
            files = [ FlagStats.flag_stats ]
    }

    # Read Names And Lengths
    call AM.ReadNamesAndLengths { input: bam = bam_file }
    call FF.FinalizeToDir as FinalizeReadNamesAndLengths {
        input:
            outdir = base_metrics_out_dir,
            files = [ ReadNamesAndLengths.read_names_and_lengths ]
    }

    output {
        String read_metrics_path = read_metrics_dir
        String depth_metrics_path = depth_metrics_dir

        # Commented out for now because mos-depth doesn't handle contigs with dashes in their names
        # https://github.com/brentp/mosdepth/issues/115
#        Array[Array[String]] chrs = MakeChrIntervalList.chrs
#
#        Array[File] full_dist      = MosDepth.full_dist
#        Array[File] global_dist    = MosDepth.global_dist
#        Array[File] region_dist    = MosDepth.region_dist
#        Array[File] regions        = MosDepth.regions
#        Array[File] regions_csi    = MosDepth.regions_csi
#        Array[File] quantized_dist = MosDepth.quantized_dist
#        Array[File] quantized      = MosDepth.quantized
#        Array[File] quantized_csi  = MosDepth.quantized_csi
#
#        Array[File] coverage_summary = SummarizeDepth.cov_summary

        File flag_stats = FlagStats.flag_stats

        File sam_stats_raw_stats = SamtoolsStats.raw_stats
        File sam_stats_summary_stats = SamtoolsStats.summary_stats
        File sam_stats_first_frag_qual_stats = SamtoolsStats.first_frag_qual
        File sam_stats_last_frag_qual_stats = SamtoolsStats.last_frag_qual
        File sam_stats_first_frag_gc_content_stats = SamtoolsStats.first_frag_gc_content
        File sam_stats_last_frag_gc_content_stats = SamtoolsStats.last_frag_gc_content
        File sam_stats_agct_content_per_cycle_stats = SamtoolsStats.acgt_content_per_cycle
        File sam_stats_insert_size_stats = SamtoolsStats.insert_size
        File sam_stats_read_length_dist_stats = SamtoolsStats.read_length_dist
        File sam_stats_indel_dist_stats = SamtoolsStats.indel_distribution
        File sam_stats_indels_per_cycle_stats = SamtoolsStats.indels_per_cycle
        File sam_stats_coverage_dist_stats = SamtoolsStats.coverage_distribution
        File sam_stats_gc_depth_stats = SamtoolsStats.gc_depth

        File read_names_and_lengths = ReadNamesAndLengths.read_names_and_lengths
    }
}
