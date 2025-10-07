version 1.0

import "../../../tasks/Utility/PBUtils.wdl" as PB
import "../../../tasks/Alignment/AlignReads.wdl" as AR
import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Preprocessing/Longbow.wdl" as Longbow
import "../../../tasks/QC/AlignedMetrics.wdl" as AM
import "../../../tasks/Visualization/NanoPlot.wdl" as NP
import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/Transcriptomics/MASSeq.wdl" as MAS

workflow PBFlowcell {

    meta {
        description: "The workflow performs the alignment of an SMRT cell's worth of data to a reference. For genomic sequencing data, the workflow also optionally performs CCS correction if the data is from a CCS library but did not get corrected on-instrument. For MAS-seq transcriptome data, this workflow will determine the most likely MAS-seq model, then it will use that model to annotate, segment, and filter the CCS reads. These CCS reads will then be aligned to the reference in trascriptome alignemnt mode. Note: Currently the MAS-seq workflow separates CLR reads, but does not process them."
    }
    parameter_meta {
        bam:                "GCS path to raw subread bam"
        ccs_report_txt:     "GCS path to CCS report txt, required if on-instrument corrected, otherwise CCS is run in this workflow for CCS libraries"
        pbi:                "GCS path to pbi index for raw subread bam"
        ref_map_file:       "table indicating reference sequence and auxillary file locations"

        SM:                 "the value to place in the BAM read group's SM field"
        LB:                 "the value to place in the BAM read group's LB (library) field"

        num_shards:         "number of shards into which fastq files should be batched"
        experiment_type:    "type of experiment run (CLR, CCS, ISOSEQ, MASSEQ)"
        dir_prefix:         "directory prefix for output files"

        mas_seq_model:      "Longbow model to use for MAS-seq data."

        DEBUG_MODE:         "[default valued] enables debugging tasks / subworkflows (default: false)"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        File bam
        File pbi
        File? ccs_report_txt

        String SM
        String LB

        File ref_map_file

        Boolean drop_per_base_N_pulse_tags = true

        Int? num_shards
        String experiment_type
        String dir_prefix

        String gcs_out_root_dir

        String? mas_seq_model

        Boolean validate_shards = false

        Boolean DEBUG_MODE = false
    }

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as WdlExecutionStartTimestamp { input: }

    Map[String, String] ref_map = read_map(ref_map_file)

    DataTypeParameters clr    = { 'num_shards': select_first([num_shards, 100]), 'map_preset': 'SUBREAD' }
    DataTypeParameters ccs    = { 'num_shards': select_first([num_shards, 100]), 'map_preset': 'CCS'     }
    DataTypeParameters isoseq = { 'num_shards': select_first([num_shards,  50]), 'map_preset': 'ISOSEQ'  }
    DataTypeParameters masseq = { 'num_shards': select_first([num_shards, 100]), 'map_preset': 'ISOSEQ'  }
    Map[String, DataTypeParameters] data_presets = { 'CLR': clr, 'CCS': ccs, 'ISOSEQ': isoseq, 'MASSEQ': masseq }

    String outdir = if DEBUG_MODE then sub(gcs_out_root_dir, "/$", "") + "/PBFlowcell/~{dir_prefix}/" + WdlExecutionStartTimestamp.timestamp_string else sub(gcs_out_root_dir, "/$", "") + "/PBFlowcell/~{dir_prefix}"

    call PB.GetRunInfo as GetRunInfo { input: bam = bam, SM = SM }
    String PU = GetRunInfo.run_info['PU']

    call Utils.GetRawReadGroup as GetRawReadGroup { input: gcs_bam_path = bam }

    if (experiment_type != "CLR" && GetRunInfo.is_corrected) {
        if (!defined(ccs_report_txt)) {
            call Utils.StopWorkflow as lack_ccs_report {input: reason = "Provided BAM is on-instrument CCS-corrected, but lacks the companion CCS report."}
        }
    }

    call Utils.ComputeAllowedLocalSSD as Guess {input: intended_gb = ceil(size(bam, "GB") * if experiment_type=='CLR' then 4 else 3 + size(pbi, "GB"))}
    call Utils.RandomZoneSpewer as arbitrary {input: num_of_zones = 3}

    # break one raw BAM into fixed number of shards
    call PB.ShardLongReads as ShardLongReads {
        input:
            unaligned_bam = bam,
            unaligned_pbi = pbi,
            num_shards = data_presets[experiment_type].num_shards,
            drop_per_base_N_pulse_tags = drop_per_base_N_pulse_tags,
            num_ssds = Guess.numb_of_local_ssd,
            zones = arbitrary.zone_string
    }

    # TODO: Peeking should only be done on CCS reads (if available)
    # for MAS-seq data, automatically detect the array model to use
    if (experiment_type == "MASSEQ" && !defined(mas_seq_model) ) {
        Int peek_size = 200
        # TODO: take first (5 * peek_size) reads from input bam, CCS correct them, feed those ccs corrected reads into peek
        call Longbow.Peek as Longbow_Peek { input: bam = ShardLongReads.unmapped_shards[0], n=peek_size }
    }

    String chosen_mas_seq_model = if (experiment_type == "MASSEQ") then select_first([mas_seq_model, Longbow_Peek.model]) else ""

    # then perform correction and alignment on each of the shard
    scatter (unmapped_shard in ShardLongReads.unmapped_shards) {
        # sometimes we see the sharded bams mising EOF marker, use this as
        if (validate_shards) {call Utils.CountBamRecords as ValidateShard {input: bam = unmapped_shard}}
        if (experiment_type != "CLR") {
            if (!GetRunInfo.is_corrected) { call PB.CCS as CCS { input: subreads = unmapped_shard } }

            File ccs_corrected_bam = select_first([CCS.consensus, unmapped_shard])

            if (experiment_type != "MASSEQ") {
                call PB.ExtractHifiReads as ExtractHifiReads {
                    input:
                        bam = ccs_corrected_bam,
                        sample_name = SM,
                        library     = LB
                }
            }
        }

        # MASSEQ Specific processing:
        if (experiment_type == "MASSEQ") {
            Float min_read_quality = 0.0

            # Split our reads into CCS and CLR:
            call Utils.Bamtools as ExtractCcsReads {
                input:
                    bamfile = select_first([ccs_corrected_bam]),
                    prefix = SM + ".ccs_reads",
                    cmd = "filter",
                    args = '-tag "rq":">=' + min_read_quality + '"',
            }
            call Utils.Bamtools as ExtractClrReads {
                input:
                    bamfile = select_first([ccs_corrected_bam]),
                    prefix = SM + ".clr_reads",
                    cmd = "filter",
                    args = '-tag "rq":"<' + min_read_quality + '"',
            }

            # Process the CCS reads:
            call Longbow.Process as LongbowProcessCCS { input: bam = ExtractCcsReads.bam, model = chosen_mas_seq_model, prefix = SM + "_" + LB }

            # Align the CCS reads:
            Array[String] masseq_tags_to_preserve =  [ "CB", "JB", "JC", "JD", "JF", "JX", "RC", "RG", "SG", "XA", "XB", "XC", "XF", "XM", "XN", "XQ", "XU", "YC", "YG", "YK", "YN", "YP", "YQ", "YS", "YV", "ZS", "ZU", "ec", "fn", "ic", "im", "is", "it", "np", "pz", "rn", "rq", "sn", "we", "ws", "zm" ]

            # Align CCS reads to the genome:
            # TODO: Debug memory requirement here:
            Int alignment_memory_gb = 32

            call AR.Minimap2 as AlignMasSeqCCSReads {
                input:
                    reads      = [ LongbowProcessCCS.extracted_bam ],
                    ref_fasta  = ref_map['fasta'],
                    RG = GetRawReadGroup.rg,
                    library = LB,
                    tags_to_preserve = masseq_tags_to_preserve,
                    map_preset = "splice:hq",
                    runtime_attr_override = object { mem_gb: alignment_memory_gb }
            }
            call Utils.BamToFastq as MasseqCCSBamToFastq { input: bam = LongbowProcessCCS.extracted_bam, prefix = basename(LongbowProcessCCS.extracted_bam, ".bam") }

            # Fix MAS-ISO-seq tags that are position-based:
            call Longbow.TagFix as FixMasSeqCcsReadTagsPostAlignment {
                input:
                    bam = AlignMasSeqCCSReads.aligned_bam,
                    prefix = SM + ".isoform_reads"
            }

            call MAS.RenameSingleCellBamTagsForMasIsoSeqV0 as RenameSingleCellBamTagsForMasIsoSeqV0 {
                input:
                    bam = FixMasSeqCcsReadTagsPostAlignment.tag_fixed_bam,
                    prefix =  SM + ".isoform_reads.renamed_tags"
            }
        }

        # All other types of libraries:
        if (experiment_type != "MASSEQ") {
            File unaligned_bam = select_first([ExtractHifiReads.hifi_bam, CCS.consensus, unmapped_shard])

            call PB.Align as AlignReads {
                input:
                    bam         = unaligned_bam,
                    ref_fasta   = ref_map['fasta'],
                    sample_name = SM,
                    library     = LB,
                    map_preset  = data_presets[experiment_type].map_preset,
                    drop_per_base_N_pulse_tags = drop_per_base_N_pulse_tags
            }

            call Utils.BamToFastq as BamToFastq { input: bam = unaligned_bam, prefix = basename(unaligned_bam, ".bam") }
        }

        # Resolve MAS-seq vs non-MAS-seq data so we can more easily reference it later:
        File unaligned_reads_bam = select_first([LongbowProcessCCS.extracted_bam, unaligned_bam])
        File aligned_reads_bam = select_first([RenameSingleCellBamTagsForMasIsoSeqV0.bam_out, AlignReads.aligned_bam])
        File reads_fastq = select_first([MasseqCCSBamToFastq.reads_fq, BamToFastq.reads_fq])
    }

    call Utils.MergeFastqs as MergeAllFastqs { input: fastqs = reads_fastq }

    # Merge the corrected per-shard BAM/report into one, corresponding to one raw input BAM
    call Utils.MergeBams as MergeAlignedReads { input: bams = aligned_reads_bam, prefix = PU }
    call PB.PBIndex as IndexAlignedReads { input: bam = MergeAlignedReads.merged_bam }

    call AM.AlignedMetrics as PerFlowcellMetrics {
        input:
            aligned_bam    = MergeAlignedReads.merged_bam,
            aligned_bai    = MergeAlignedReads.merged_bai,
            ref_fasta      = ref_map['fasta'],
            ref_dict       = ref_map['dict'],
            gcs_output_dir = outdir + "/metrics"
    }

    # For finalization, we use a keyfile that should be generated last:
    File keyfile = PerFlowcellMetrics.raw_chr_intervals

    # Merge corrected, unaligned reads:
    String cdir = outdir + "/reads/ccs/unaligned"
    if (experiment_type != "CLR") {
        call Utils.MergeBams as MergeCCSUnalignedReads { input: bams = unaligned_reads_bam, prefix = "~{PU}.reads" }
        call PB.PBIndex as IndexCCSUnalignedReads { input: bam = MergeCCSUnalignedReads.merged_bam }

        call FF.FinalizeToFile as FinalizeCCSUnalignedBam { input: outdir = cdir, file = MergeCCSUnalignedReads.merged_bam, keyfile = keyfile }
        call FF.FinalizeToFile as FinalizeCCSUnalignedBai { input: outdir = cdir, file = MergeCCSUnalignedReads.merged_bai, keyfile = keyfile }
        call FF.FinalizeToFile as FinalizeCCSUnalignedPbi {
            input:
                outdir = cdir,
                file = IndexCCSUnalignedReads.pbi,
                name = basename(MergeCCSUnalignedReads.merged_bam) + ".pbi",
                keyfile = keyfile
        }
    }

    # Merge CCS Reports:
    if (experiment_type != "CLR" && !GetRunInfo.is_corrected) {
        call PB.MergeCCSReports as MergeCCSReports { input: reports = select_all(CCS.report), prefix = PU }
        call PB.SummarizeCCSReport as SummarizeCCSReport { input: report = MergeCCSReports.report }

        call FF.FinalizeToFile as FinalizeCCSReport { input: outdir = cdir, file = MergeCCSReports.report, keyfile = keyfile }
    }

    if (experiment_type == 'MASSEQ') {
        String longbow_stats_dir = outdir + "/metrics/longbow"

        # Merge all longbow intermediate files:
        call Utils.MergeBams as MergeLongbowAnnotatedCCSBams { input: bams = select_all(LongbowProcessCCS.annotated_bam), prefix = "~{PU}.reads.longbow_annotated" }
        call Utils.MergeBams as MergeLongbowFilteredCCSBams { input: bams = select_all(LongbowProcessCCS.filtered_bam), prefix = "~{PU}.reads.longbow_annotated_filter_passed" }

        # Overall CCS Stats:
        call Longbow.Stats as LongbowOverallCCSStats { input: bam = MergeLongbowAnnotatedCCSBams.merged_bam, prefix="overall_ccs" }
        call FF.FinalizeToDir as FinalizeOverallCCSStatsPNGs { input: outdir = longbow_stats_dir, files = LongbowOverallCCSStats.pngs, keyfile = keyfile }
        call FF.FinalizeToDir as FinalizeOverallCCSStatsSVGs { input: outdir = longbow_stats_dir, files = LongbowOverallCCSStats.svgs, keyfile = keyfile }
        call FF.FinalizeToFile as FinalizeOverallCCSStatsSummary { input: outdir = longbow_stats_dir, file = LongbowOverallCCSStats.summary, keyfile = keyfile }

        # Filtered CCS Stats:
        call Longbow.Stats as LongbowFilteredCCSStats { input: bam = MergeLongbowFilteredCCSBams.merged_bam, prefix="filtered_ccs" }
        call FF.FinalizeToDir as FinalizeFilteredCCSStatsPNGs { input: outdir = longbow_stats_dir, files = LongbowFilteredCCSStats.pngs, keyfile = keyfile }
        call FF.FinalizeToDir as FinalizeFilteredCCSStatsSVGs { input: outdir = longbow_stats_dir, files = LongbowFilteredCCSStats.svgs, keyfile = keyfile }
        call FF.FinalizeToFile as FinalizeFilteredCCSStatsSummary { input: outdir = longbow_stats_dir, file = LongbowFilteredCCSStats.summary, keyfile = keyfile }

        # Some debugging outputs:
        String longbow_intermediate_reads_dir = outdir + "/reads/longbow/intermediate"

        if (DEBUG_MODE) {
            call Utils.MergeBams as MergeLongbowSegmentedCCSBams { input: bams = select_all(LongbowProcessCCS.segmented_bam), prefix = "~{PU}.reads.longbow_annotated_segmented" }
            call Utils.MergeBams as MergeLongbowFilterFailedCCSBams { input: bams = select_all(LongbowProcessCCS.filter_failed_bam), prefix = "~{PU}.reads.longbow_annotated_filter_failed" }
            call Utils.MergeBams as MergeLongbowExtractedCCSBams { input: bams = select_all(LongbowProcessCCS.extracted_bam), prefix = "~{PU}.reads.longbow_annotated_segmented_filter_passed_corrected_extracted" }

            call FF.FinalizeToFile as FinalizeLongbowAnnotatedBam { input: outdir = longbow_intermediate_reads_dir, file = MergeLongbowAnnotatedCCSBams.merged_bam, keyfile = keyfile }
            call FF.FinalizeToFile as FinalizeLongbowSegmentedBam { input: outdir = longbow_intermediate_reads_dir, file = MergeLongbowSegmentedCCSBams.merged_bam, keyfile = keyfile }
            call FF.FinalizeToFile as FinalizeLongbowFilteredBam { input: outdir = longbow_intermediate_reads_dir, file = MergeLongbowFilteredCCSBams.merged_bam, keyfile = keyfile }
            call FF.FinalizeToFile as FinalizeLongbowFilterFailedBam { input: outdir = longbow_intermediate_reads_dir, file = MergeLongbowFilterFailedCCSBams.merged_bam, keyfile = keyfile }

            call FF.FinalizeToFile as FinalizeLongbowExtractedBam { input: outdir = longbow_intermediate_reads_dir, file = MergeLongbowExtractedCCSBams.merged_bam, keyfile = keyfile }
            call FF.FinalizeToDir as FinalizeLongbowIntermediateBamIndices { input: outdir = longbow_intermediate_reads_dir, files = [ MergeLongbowAnnotatedCCSBams.merged_bai, MergeLongbowSegmentedCCSBams.merged_bai, MergeLongbowFilteredCCSBams.merged_bai, MergeLongbowFilterFailedCCSBams.merged_bai, MergeLongbowExtractedCCSBams.merged_bai ], keyfile = keyfile}

            # Only merge and finalize these files if they were created:
            if (chosen_mas_seq_model != "mas_15_bulk_10x5p_single_internal" && chosen_mas_seq_model != "mas_15_bulk_10x5p_single_internal") {
                call Utils.MergeBams as MergeLongbowCorrectedCCSBams { input: bams = select_all(LongbowProcessCCS.corrected_bam), prefix = "~{PU}.reads.longbow_annotated_segmented_filter_passed_corrected" }
                call Utils.MergeBams as MergeLongbowUncorrectableCCSBams { input: bams = select_all(LongbowProcessCCS.uncorrectable_bam), prefix = "~{PU}.reads.longbow_annotated_segmented_filter_passed_uncorrectable" }
                call FF.FinalizeToFile as FinalizeLongbowCorrectedBam { input: outdir = longbow_intermediate_reads_dir, file = MergeLongbowCorrectedCCSBams.merged_bam, keyfile = keyfile }
                call FF.FinalizeToFile as FinalizeLongbowUncorrectedBam { input: outdir = longbow_intermediate_reads_dir, file = MergeLongbowUncorrectableCCSBams.merged_bam, keyfile = keyfile }
                call FF.FinalizeToDir as FinalizeLongbowIntermediateSingleCellBamIndices { input: outdir = longbow_intermediate_reads_dir, files = [ MergeLongbowCorrectedCCSBams.merged_bai, MergeLongbowUncorrectableCCSBams.merged_bai], keyfile = keyfile}
            }
        }

        # Merge and finalize CLR reads:
        call Utils.MergeBams as MergeMASSeqCLRReads { input: bams = select_all(ExtractClrReads.bam), prefix = "~{PU}.CLR_reads" }
        String longbow_intermediate_clr_reads_dir = outdir + "/reads/clr"
        call FF.FinalizeToFile as FinalizeMasSeqUnprocessedCLRReads { input: outdir = longbow_intermediate_clr_reads_dir, file = MergeMASSeqCLRReads.merged_bam, keyfile = keyfile }
        call FF.FinalizeToFile as FinalizeMasSeqUnprocessedCLRReadsBai { input: outdir = longbow_intermediate_clr_reads_dir, file = MergeMASSeqCLRReads.merged_bai, keyfile = keyfile }

        # Finalize MAS-seq logs / stats:
        String longbow_logs_dir = outdir + "/logs"
        scatter (i_1 in range(length(select_all(LongbowProcessCCS.correct_log)))) {

            # Only finalize these files if they were created:
            if (chosen_mas_seq_model != "mas_15_bulk_10x5p_single_internal" && chosen_mas_seq_model != "mas_15_bulk_10x5p_single_internal") {
                call FF.FinalizeToFile as FinalizeMasSeqLongbowCorrectStats { input: outdir = longbow_stats_dir, file = select_first([LongbowProcessCCS.correct_stats[i_1]]), name = "longbow_correct_stats.shard_" + i_1 + ".txt", keyfile = keyfile }
            }

            # Some more debugging outputs:
            if (DEBUG_MODE) {
                call FF.FinalizeToFile as FinalizeLongbowUmiAdjustmentLog { input: outdir = longbow_logs_dir + "/longbow_umi_adjustment", file = select_first([LongbowProcessCCS.umi_adjustment_log[i_1]]), name = "longbow_umi_adjustment.shard_" + i_1 + ".txt", keyfile = keyfile }
                # Only finalize these files if they were created:
                if (chosen_mas_seq_model != "mas_15_bulk_10x5p_single_internal" && chosen_mas_seq_model != "mas_15_bulk_10x5p_single_internal") {
                    File correct_logs_shard = select_first([LongbowProcessCCS.correct_log[i_1]])
                    call FF.FinalizeToFile as FinalizeLongbowCorrectLog { input: outdir = longbow_logs_dir + "/longbow_correct", file = correct_logs_shard, name = "longbow_correct_stats.shard_" + i_1 + ".txt", keyfile = keyfile }
                }
            }
        }

        # TODO: Run a Jupyter notebook to create a report for us in an HTML file:
    }

    # Get the correct PBI file on which to calculate stats:
    File subreads_pbi_file = if (experiment_type == "MASSEQ") then select_first([IndexCCSUnalignedReads.pbi]) else pbi
    call PB.SummarizePBI as SummarizeSubreadsPBI   { input: pbi = subreads_pbi_file, runtime_attr_override = { 'mem_gb': 72 } }

    # Collect stats on aligned reads:
    call PB.SummarizePBI as SummarizeAlignedQ5PBI  { input: pbi = IndexAlignedReads.pbi, qual_threshold = 5 }
    call PB.SummarizePBI as SummarizeAlignedQ7PBI  { input: pbi = IndexAlignedReads.pbi, qual_threshold = 7 }
    call PB.SummarizePBI as SummarizeAlignedQ10PBI { input: pbi = IndexAlignedReads.pbi, qual_threshold = 10 }
    call PB.SummarizePBI as SummarizeAlignedQ12PBI { input: pbi = IndexAlignedReads.pbi, qual_threshold = 12 }
    call PB.SummarizePBI as SummarizeAlignedQ15PBI { input: pbi = IndexAlignedReads.pbi, qual_threshold = 15 }

    call NP.NanoPlotFromBam as NanoPlotFromBam { input: bam = MergeAlignedReads.merged_bam, bai = MergeAlignedReads.merged_bai }
    call Utils.ComputeGenomeLength as ComputeGenomeLength { input: fasta = ref_map['fasta'] }

    # Finalize data
    String dir = outdir + "/" + if (experiment_type != "CLR") then "reads/ccs/aligned" else "reads/subreads/aligned"

    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = dir, file = MergeAlignedReads.merged_bam, keyfile = keyfile }
    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = dir, file = MergeAlignedReads.merged_bai, keyfile = keyfile }
    call FF.FinalizeToFile as FinalizeAlignedPbi {
        input:
            outdir = dir,
            file = IndexAlignedReads.pbi,
            name = basename(MergeAlignedReads.merged_bam) + ".pbi",
            keyfile = keyfile
    }

    String fqdir = outdir + "/" + if (experiment_type != "CLR") then "reads/ccs/unaligned" else "reads/subreads/unaligned"
    call FF.FinalizeToFile as FinalizeFastq {
        input:
            outdir = fqdir,
            file = MergeAllFastqs.merged_fastq,
            name = basename(MergeAlignedReads.merged_bam, ".bam") + ".fq.gz",
            keyfile = keyfile
    }

    # Set our library type for our outputs:
    String lib_type = if (experiment_type == "MASSEQ") then chosen_mas_seq_model else experiment_type

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good:
    call FF.WriteCompletionFile as WriteCompletionFile {
        input:
            outdir = outdir + "/",
            keyfile = keyfile
    }

    output {
        # Flowcell stats
        File? ccs_report = FinalizeCCSReport.gcs_path
        Float? ccs_zmws_input = SummarizeCCSReport.zmws_input
        Float? ccs_zmws_pass_filters = SummarizeCCSReport.zmws_pass_filters
        Float? ccs_zmws_fail_filters = SummarizeCCSReport.zmws_fail_filters
        Float? ccs_zmws_shortcut_filters = SummarizeCCSReport.zmws_shortcut_filters
        Float? ccs_zmws_pass_filters_pct = SummarizeCCSReport.zmws_pass_filters_pct
        Float? ccs_zmws_fail_filters_pct = SummarizeCCSReport.zmws_fail_filters_pct
        Float? ccs_zmws_shortcut_filters_pct = SummarizeCCSReport.zmws_shortcut_filters_pct

        Float polymerase_read_length_mean = SummarizeSubreadsPBI.results['polymerase_mean']
        Float polymerase_read_length_N50 = SummarizeSubreadsPBI.results['polymerase_n50']

        Float subread_read_length_mean = SummarizeSubreadsPBI.results['subread_mean']
        Float subread_read_length_N50 = SummarizeSubreadsPBI.results['subread_n50']

        # Unaligned reads
        File? fq = FinalizeFastq.gcs_path

        # Unaligned BAM file
        File? ccs_bam = FinalizeCCSUnalignedBam.gcs_path
        File? ccs_pbi = FinalizeCCSUnalignedPbi.gcs_path

        # Aligned BAM file
        File aligned_bam = FinalizeAlignedBam.gcs_path
        File aligned_bai = FinalizeAlignedBai.gcs_path
        File aligned_pbi = FinalizeAlignedPbi.gcs_path

        # Unaligned read stats
        Float num_reads = SummarizeSubreadsPBI.results['reads']
        Float num_bases = SummarizeSubreadsPBI.results['bases']
        Float raw_est_fold_cov = SummarizeSubreadsPBI.results['bases']/ComputeGenomeLength.length

        Float read_length_mean = SummarizeSubreadsPBI.results['subread_mean']
        Float read_length_median = SummarizeSubreadsPBI.results['subread_median']
        Float read_length_stdev = SummarizeSubreadsPBI.results['subread_stdev']
        Float read_length_N50 = SummarizeSubreadsPBI.results['subread_n50']

        Float read_qual_mean = SummarizeSubreadsPBI.results['mean_qual']
        Float read_qual_median = SummarizeSubreadsPBI.results['median_qual']

        Float num_reads_Q5 = SummarizeAlignedQ5PBI.results['reads']
        Float num_reads_Q7 = SummarizeAlignedQ7PBI.results['reads']
        Float num_reads_Q10 = SummarizeAlignedQ10PBI.results['reads']
        Float num_reads_Q12 = SummarizeAlignedQ12PBI.results['reads']
        Float num_reads_Q15 = SummarizeAlignedQ15PBI.results['reads']

        # Aligned read stats
        Float aligned_num_reads = NanoPlotFromBam.stats_map['number_of_reads']
        Float aligned_num_bases = NanoPlotFromBam.stats_map['number_of_bases_aligned']
        Float aligned_frac_bases = NanoPlotFromBam.stats_map['fraction_bases_aligned']
        Float aligned_est_fold_cov = NanoPlotFromBam.stats_map['number_of_bases_aligned']/ComputeGenomeLength.length

        Float aligned_read_length_mean = NanoPlotFromBam.stats_map['mean_read_length']
        Float aligned_read_length_median = NanoPlotFromBam.stats_map['median_read_length']
        Float aligned_read_length_stdev = NanoPlotFromBam.stats_map['read_length_stdev']
        Float aligned_read_length_N50 = NanoPlotFromBam.stats_map['n50']

        Float average_identity = NanoPlotFromBam.stats_map['average_identity']
        Float median_identity = NanoPlotFromBam.stats_map['median_identity']

        String library_type = lib_type

        # MAS-seq outputs:
        # Overall CCS Stats:
        String? longbow_overall_ccs_stats_plots_png = FinalizeOverallCCSStatsPNGs.gcs_dir
        String? longbow_overall_ccs_stats_plots_svg = FinalizeOverallCCSStatsSVGs.gcs_dir
        File? longbow_overall_ccs_stats = FinalizeOverallCCSStatsSummary.gcs_path

        # Filtered CCS Stats:
        String? longbow_filtered_ccs_stats_plots_png = FinalizeFilteredCCSStatsPNGs.gcs_dir
        String? longbow_filtered_ccs_stats_plots_svg = FinalizeFilteredCCSStatsSVGs.gcs_dir
        File? longbow_filtered_ccs_stats = FinalizeFilteredCCSStatsSummary.gcs_path
    }
}
