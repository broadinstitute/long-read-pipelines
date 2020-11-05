version 1.0

##########################################################################################
## A workflow that performs CCS correction and variant calling on PacBio HiFi reads from a
## single flow cell. The workflow shards the subreads into clusters and performs CCS in
## parallel on each cluster.  Error-corrected reads are then variant-called.  A number of
## metrics and figures are produced along the way.
##########################################################################################

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/CallSVs.wdl" as SV
import "tasks/Figures.wdl" as FIG
import "tasks/Finalize.wdl" as FF
import "tasks/CallSmallVariants.wdl" as SMV

workflow PBCCSWholeGenome {
    input {
        Array[File] bams
        File ref_map_file

        String participant_name
        Int num_shards = 300
        Boolean extract_uncorrected_reads = false

        String? gcs_out_root_dir
    }

    parameter_meta {
        bams:                      "GCS path to raw subreads or CCS data"
        ref_map_file:              "Table indicating reference sequence and auxillary file locations"

        participant_name:          "name of the participant from whom these samples were obtained"
        num_shards:                "[default-valued] number of sharded BAMs to create (tune for performance)"
        extract_uncorrected_reads: "[default-valued] extract reads that were not CCS-corrected to a separate file"

        gcs_out_root_dir:          "[optional] GCS bucket to store the corrected/uncorrected reads, variants, and metrics files"
    }

    call Utils.GetDefaultDir

    String outdir = sub(select_first([gcs_out_root_dir, GetDefaultDir.path]), "/$", "") + "/" + participant_name

    Map[String, String] ref_map = read_map(ref_map_file)

    # scatter over all sample BAMs
    scatter (bam in bams) {
        File pbi = sub(bam, ".bam$", ".bam.pbi")

        call PB.GetRunInfo { input: bam = bam }
        String ID = GetRunInfo.run_info["PU"]

        # break one raw BAM into fixed number of shards
        call PB.ShardLongReads { input: unaligned_bam = bam, unaligned_pbi = pbi, num_shards = num_shards }

        # then perform correction and alignment on each of the shard
        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads }

            if (extract_uncorrected_reads) {
                call PB.ExtractUncorrectedReads { input: subreads = subreads, consensus = CCS.consensus }

                call PB.Align as AlignUncorrected {
                    input:
                        bam         = ExtractUncorrectedReads.uncorrected,
                        ref_fasta   = ref_map['fasta'],
                        sample_name = participant_name,
                        map_preset  = "SUBREAD",
                        runtime_attr_override = { 'mem_gb': 64 }
                }
            }

            call PB.Align as AlignCorrected {
                input:
                    bam         = CCS.consensus,
                    ref_fasta   = ref_map['fasta'],
                    sample_name = participant_name,
                    map_preset  = "CCS"
            }
        }

        # merge the corrected per-shard BAM/report into one, corresponding to one raw input BAM
        call Utils.MergeBams as MergeCorrected { input: bams = AlignCorrected.aligned_bam, prefix = "~{participant_name}.~{ID}.corrected" }

        if (length(select_all(AlignUncorrected.aligned_bam)) > 0) {
            call Utils.MergeBams as MergeUncorrected {
                input:
                    bams = select_all(AlignUncorrected.aligned_bam),
                    prefix = "~{participant_name}.~{ID}.uncorrected"
            }
        }

        # compute alignment metrics
        call AM.AlignedMetrics as PerFlowcellMetrics {
            input:
                aligned_bam    = MergeCorrected.merged_bam,
                aligned_bai    = MergeCorrected.merged_bai,
                ref_fasta      = ref_map['fasta'],
                ref_dict       = ref_map['dict'],
                ref_flat       = ref_map['flat'],
                dbsnp_vcf      = ref_map['dbsnp_vcf'],
                dbsnp_tbi      = ref_map['dbsnp_tbi'],
                metrics_locus  = ref_map['metrics_locus'],
                gcs_output_dir = outdir + "/metrics/per_flowcell/" + ID
        }

        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report, prefix = "~{participant_name}.~{ID}" }

        call FF.FinalizeToDir as FinalizeCCSReport {
            input:
                files = [ MergeCCSReports.report ],
                outdir = outdir + "/metrics/per_flowcell/" + ID + "/ccs_metrics"
        }
    }

    # gather across (potential multiple) input raw BAMs
    if (length(bams) > 1) {
        call Utils.MergeBams as MergeAllCorrected { input: bams = MergeCorrected.merged_bam, prefix = "~{participant_name}.corrected" }
        call PB.MergeCCSReports as MergeAllCCSReports { input: reports = MergeCCSReports.report }

        if (length(select_all(MergeUncorrected.merged_bam)) > 0) {
            call Utils.MergeBams as MergeAllUncorrected {
                input:
                    bams = select_all(MergeUncorrected.merged_bam),
                    prefix = "~{participant_name}.uncorrected"
            }
        }
    }

    File ccs_bam = select_first([ MergeAllCorrected.merged_bam, MergeCorrected.merged_bam[0] ])
    File ccs_bai = select_first([ MergeAllCorrected.merged_bai, MergeCorrected.merged_bai[0] ])
    File ccs_report = select_first([ MergeAllCCSReports.report, MergeCCSReports.report[0] ])

    if (extract_uncorrected_reads) {
        File? uncorrected_bam = select_first([ MergeAllUncorrected.merged_bam, MergeUncorrected.merged_bam[0] ])
        File? uncorrected_bai = select_first([ MergeAllUncorrected.merged_bai, MergeUncorrected.merged_bai[0] ])
    }

    # compute alignment metrics
    call AM.AlignedMetrics as PerSampleMetrics {
        input:
            aligned_bam    = ccs_bam,
            aligned_bai    = ccs_bai,
            ref_fasta      = ref_map['fasta'],
            ref_dict       = ref_map['dict'],
            ref_flat       = ref_map['flat'],
            dbsnp_vcf      = ref_map['dbsnp_vcf'],
            dbsnp_tbi      = ref_map['dbsnp_tbi'],
            metrics_locus  = ref_map['metrics_locus'],
            gcs_output_dir = outdir + "/metrics/combined/" + participant_name
    }

    # call SVs
    call SV.CallSVs as CallSVs {
        input:
            bam               = ccs_bam,
            bai               = ccs_bai,

            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
            tandem_repeat_bed = ref_map['tandem_repeat_bed'],

            preset            = "hifi"
    }

    # call SNVs and small indels
    call SMV.CallSmallVariants as CallSmallVariants {
        input:
            bam               = ccs_bam,
            bai               = ccs_bai,

            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
            ref_dict          = ref_map['dict'],
    }

    ##########
    # store the results into designated bucket
    ##########

    call FF.FinalizeToDir as FinalizeSVs {
        input:
            files = [ CallSVs.pbsv_vcf, CallSVs.sniffles_vcf, CallSVs.svim_vcf, CallSVs.cutesv_vcf ],
            outdir = outdir + "/variants"
    }

    call FF.FinalizeToDir as FinalizeSmallVariants {
        input:
            files = [ CallSmallVariants.longshot_vcf, CallSmallVariants.longshot_tbi ],
            outdir = outdir + "/variants"
    }

    call FF.FinalizeToDir as FinalizeMergedRuns {
        input:
            files = [ ccs_bam, ccs_bai ],
            outdir = outdir + "/alignments"
    }

    call FF.FinalizeToDir as FinalizeMergedCCSReport {
        input:
            files = [ ccs_report ],
            outdir = outdir + "/metrics/combined/" + participant_name + "/ccs_metrics"
    }

    if (extract_uncorrected_reads) {
        call FF.FinalizeToDir as FinalizeMergedUncorrectedRuns {
            input:
                files = select_all([ uncorrected_bam, uncorrected_bai ]),
                outdir = outdir + "/alignments"
        }
    }

    output {
        # BAMs
        File corrected_bam = ccs_bam
        File corrected_bai = ccs_bai

        # SVs
        File pbsv_vcf = CallSVs.pbsv_vcf
        File sniffles_vcf = CallSVs.sniffles_vcf
        File svim_vcf = CallSVs.svim_vcf
        File cutesv_vcf = CallSVs.cutesv_vcf

        # SNPs/indels
        File longshot_vcf = CallSmallVariants.longshot_vcf
        File longshot_tbi = CallSmallVariants.longshot_tbi

        # Per-sample metrics
        File per_sample_aligned_flag_stats = PerSampleMetrics.aligned_flag_stats

        Array[File] per_sample_coverage_full_dist      = PerSampleMetrics.coverage_full_dist
        Array[File] per_sample_coverage_global_dist    = PerSampleMetrics.coverage_global_dist
        Array[File] per_sample_coverage_region_dist    = PerSampleMetrics.coverage_region_dist
        Array[File] per_sample_coverage_regions        = PerSampleMetrics.coverage_regions
        Array[File] per_sample_coverage_regions_csi    = PerSampleMetrics.coverage_regions_csi
        Array[File] per_sample_coverage_quantized_dist = PerSampleMetrics.coverage_quantized_dist
        Array[File] per_sample_coverage_quantized      = PerSampleMetrics.coverage_quantized
        Array[File] per_sample_coverage_quantized_csi  = PerSampleMetrics.coverage_quantized_csi

        File per_sample_aligned_np_hist        = PerSampleMetrics.aligned_np_hist
        File per_sample_aligned_range_gap_hist = PerSampleMetrics.aligned_range_gap_hist
        File per_sample_aligned_zmw_hist       = PerSampleMetrics.aligned_zmw_hist
        File per_sample_aligned_prl_counts     = PerSampleMetrics.aligned_prl_counts
        File per_sample_aligned_prl_hist       = PerSampleMetrics.aligned_prl_hist
        File per_sample_aligned_prl_nx         = PerSampleMetrics.aligned_prl_nx
        File per_sample_aligned_prl_yield_hist = PerSampleMetrics.aligned_prl_yield_hist
        File per_sample_aligned_rl_counts      = PerSampleMetrics.aligned_rl_counts
        File per_sample_aligned_rl_hist        = PerSampleMetrics.aligned_rl_hist
        File per_sample_aligned_rl_nx          = PerSampleMetrics.aligned_rl_nx
        File per_sample_aligned_rl_yield_hist  = PerSampleMetrics.aligned_rl_yield_hist

        # Per-flowcell metrics
        Array[File] per_flowcell_aligned_flag_stats = PerFlowcellMetrics.aligned_flag_stats

        Array[Array[File]] per_flowcell_coverage_full_dist      = PerFlowcellMetrics.coverage_full_dist
        Array[Array[File]] per_flowcell_coverage_global_dist    = PerFlowcellMetrics.coverage_global_dist
        Array[Array[File]] per_flowcell_coverage_region_dist    = PerFlowcellMetrics.coverage_region_dist
        Array[Array[File]] per_flowcell_coverage_regions        = PerFlowcellMetrics.coverage_regions
        Array[Array[File]] per_flowcell_coverage_regions_csi    = PerFlowcellMetrics.coverage_regions_csi
        Array[Array[File]] per_flowcell_coverage_quantized_dist = PerFlowcellMetrics.coverage_quantized_dist
        Array[Array[File]] per_flowcell_coverage_quantized      = PerFlowcellMetrics.coverage_quantized
        Array[Array[File]] per_flowcell_coverage_quantized_csi  = PerFlowcellMetrics.coverage_quantized_csi

        Array[File] per_flowcell_aligned_np_hist        = PerFlowcellMetrics.aligned_np_hist
        Array[File] per_flowcell_aligned_range_gap_hist = PerFlowcellMetrics.aligned_range_gap_hist
        Array[File] per_flowcell_aligned_zmw_hist       = PerFlowcellMetrics.aligned_zmw_hist
        Array[File] per_flowcell_aligned_prl_counts     = PerFlowcellMetrics.aligned_prl_counts
        Array[File] per_flowcell_aligned_prl_hist       = PerFlowcellMetrics.aligned_prl_hist
        Array[File] per_flowcell_aligned_prl_nx         = PerFlowcellMetrics.aligned_prl_nx
        Array[File] per_flowcell_aligned_prl_yield_hist = PerFlowcellMetrics.aligned_prl_yield_hist
        Array[File] per_flowcell_aligned_rl_counts      = PerFlowcellMetrics.aligned_rl_counts
        Array[File] per_flowcell_aligned_rl_hist        = PerFlowcellMetrics.aligned_rl_hist
        Array[File] per_flowcell_aligned_rl_nx          = PerFlowcellMetrics.aligned_rl_nx
        Array[File] per_flowcell_aligned_rl_yield_hist  = PerFlowcellMetrics.aligned_rl_yield_hist
    }
}
