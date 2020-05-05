version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF
import "tasks/AlignReads.wdl" as AR
import "tasks/Cartographer.wdl" as CART
import "tasks/TranscriptAnalysis/Flair_Tasks.wdl" as ISO
import "tasks/ReadsMetrics.wdl" as RM
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/Ten_X_Tool.wdl" as TENX
import "tasks/JupyterNotebooks.wdl" as JUPYTER

workflow PB10xMasSeqSingleFlowcell {

    meta {
        description : "This workflow is designed to process data from the MASSeq protocol and produce aligned reads that are ready for downstream analysis (e.g. transcript isoform identification).  It takes in a raw PacBio run folder location on GCS and produces a folder containing the aligned reads and other processed data."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String gcs_input_dir
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/PB10xMasSeqSingleFlowcell"

        File segments_fasta = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.1/cDNA_array_8x.unique_seqs_for_cartographer.fasta"
        File boundaries_file = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.1/bounds_file_for_extraction.txt"

        File head_adapter_fasta = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.1/10x_adapter.fasta"
        File tail_adapter_fasta = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.1/tso_adapter.fasta"

        File ref_fasta =  "gs://broad-dsde-methods-long-reads/resources/references/grch38/Homo_sapiens_assembly38.fasta"
        File ref_fasta_index = "gs://broad-dsde-methods-long-reads/resources/references/grch38/Homo_sapiens_assembly38.fasta.fai"
        File ref_fasta_dict = "gs://broad-dsde-methods-long-reads/resources/references/grch38/Homo_sapiens_assembly38.dict"

        File ref_flat_annotations = "gs://broad-dsde-methods-long-reads/resources/references/grch38/refFlat.txt"

        File jupyter_template_interactive = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.1/MAS-seq_QC_report_template-interactive.ipynb"
        File jupyter_template_static = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.1/MAS-seq_QC_report_template-static.ipynb"
        File workflow_dot_file = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.1/PB10xMasSeqArraySingleFlowcell.dot"

        String? sample_name
    }

    parameter_meta {
        gcs_input_dir : "Input folder on GCS in which to search for BAM files to process."
        gcs_out_root_dir : "Root output GCS folder in which to place results of this workflow."

        segments_fasta : "FASTA file containing unique segments for which to search in the given BAM files.   These segments are used as delimiters in the reads.  Read splitting uses these delimiters and the boundaries file."
        boundaries_file : "Text file containing two comma-separated segment names from the segments_fasta on each line.  These entries define delimited sections to be extracted from the reads and treated as individual array elements."

        head_adapter_fasta : "FASTA file containing the sequence that each transcript should start with.  Typically this will be the 10x adapter sequence from the 10x library prep."
        tail_adapter_fasta : "FASTA file containing the sequence that each transcript should end with.  Typically this will be the Template Switch Oligo (TSO) sequence from the 10x library prep."

        ref_fasta : "FASTA file containing the reference sequence to which the input data should be aligned."
        ref_fasta_index : "FASTA index file for the given ref_fasta file."
        ref_fasta_dict : "Sequence dictionary file for the given ref_fasta file."

        ref_flat_annotations : "RefFlat file containing genomic annotations (used for RNASeq metrics).  This file must match the reference to which the input data are aligned."

        jupyter_template_interactive : "Jupyter notebook / ipynb file containing a template for the QC report which will contain interactive plots (such as those created with Bokeh and Plot.ly).  This should contain the same information as the jupyter_template_static file, but with static images."
        jupyter_template_static : "Jupyter notebook / ipynb file containing a template for the QC report which will contain static plots.  This should contain the same information as the jupyter_template_interactive file, but with static images."
        workflow_dot_file : "DOT file containing the representation of this WDL to be included in the QC reports.  This can be generated with womtool."

        sample_name : "[optional] The name of the sample to associate with the data in this workflow."
    }

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as WdlExecutionStartTimestamp { input: }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = gcs_input_dir }

    scatter (subread_bam in FindBams.subread_bams) {
        call PB.GetRunInfo { input: subread_bam = subread_bam }

        String SM  = select_first([sample_name, GetRunInfo.run_info["SM"]])
        String PL  = "PACBIO"
        String PU  = GetRunInfo.run_info["PU"]
        String DT  = GetRunInfo.run_info["DT"]
        String ID  = PU
        String DS  = GetRunInfo.run_info["DS"]
        String DIR = SM + "." + ID

        String RG_subreads  = "@RG\\tID:~{ID}.subreads\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
        String RG_consensus = "@RG\\tID:~{ID}.consensus\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
        String RG_array_elements = "@RG\\tID:~{ID}.array_elements\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        call Utils.ShardLongReadsWithCopy { input: unmapped_files = [ subread_bam ], num_reads_per_split = 2000000 }

        scatter (subreads in ShardLongReadsWithCopy.unmapped_shards) {

            # Call CCS on the subreads from the sequencer:
            call PB.CCS { input: subreads = subreads, preemptible_attempts = 0 }

            # Shard these reads even wider so we can make sure we don't run out of memory:
            call Utils.ShardLongReadsWithCopy as ShardCorrectedReads { input: unmapped_files = [ CCS.consensus ], num_reads_per_split = 20000 }

            scatter (corrected_shard in ShardCorrectedReads.unmapped_shards) {

                call CART.ExtractBoundedReadSectionsTask {
                    input:
                        reads_file          = corrected_shard,
                        segments_fasta      = segments_fasta,
                        boundaries_file     = boundaries_file,
                        max_read_length     = 50000,
                        mem_gb              = 16
                }
            }

            # Merge all outputs of ExtractBoundedReadSectionsTask:
            call Utils.MergeFiles as MergeArrayElementSubShards_1 { input: files_to_merge = ExtractBoundedReadSectionsTask.extracted_reads, merged_file_name = "EBR_extracted_reads.fasta" }
            call Utils.MergeFiles as MergeArrayElementRejectedReads_1 { input: files_to_merge = ExtractBoundedReadSectionsTask.rejected_reads, merged_file_name = "EBR_rejected_reads.txt" }
            call Utils.MergeFiles as MergeArrayElementMarkerAlignments_1 { input: files_to_merge = ExtractBoundedReadSectionsTask.raw_marker_alignments, merged_file_name = "EBR_marker_alignments.txt" }
            call Utils.MergeFiles as MergeArrayElementInitialSections_1 { input: files_to_merge = ExtractBoundedReadSectionsTask.initial_section_alignments, merged_file_name = "EBR_initial_sections.txt" }
            call Utils.MergeFiles as MergeArrayElementFinalSections_1 { input: files_to_merge = ExtractBoundedReadSectionsTask.final_section_alignments, merged_file_name = "EBR_final_sections.txt" }

            call AR.Minimap2 as AlignArrayElements {
                input:
                    reads      = [ MergeArrayElementSubShards_1.merged_file ],
                    ref_fasta  = ref_fasta,
                    RG         = RG_array_elements,
                    map_preset = "splice"
            }

            call AR.Minimap2 as AlignCCSReads {
                input:
                    reads      = [ CCS.consensus ],
                    ref_fasta  = ref_fasta,
                    RG         = RG_consensus,
                    map_preset = "splice"
            }

            call TENX.AnnotateBarcodesAndUMIs as AnnotateArrayElements {
                input:
                    bam_file = AlignArrayElements.aligned_bam,
                    bam_index = AlignArrayElements.aligned_bai,
                    head_adapter_fasta = head_adapter_fasta,
                    tail_adapter_fasta = tail_adapter_fasta,
                    read_end_length = 200,
                    poly_t_length = 31,
                    barcode_length = 16,
                    umi_length = 12
            }
        }

        # Merge all sharded merged outputs of ExtractBoundedReadSectionsTask:
        call Utils.MergeFiles as MergeArrayElementSubShards_2 { input: files_to_merge = MergeArrayElementSubShards_1.merged_file, merged_file_name = "EBR_extracted_reads.fasta" }
        call Utils.MergeFiles as MergeArrayElementRejectedReads_2 { input: files_to_merge = MergeArrayElementRejectedReads_1.merged_file, merged_file_name = "EBR_rejected_reads.txt" }
        call Utils.MergeFiles as MergeArrayElementMarkerAlignments_2 { input: files_to_merge = MergeArrayElementMarkerAlignments_1.merged_file, merged_file_name = "EBR_marker_alignments.txt" }
        call Utils.MergeFiles as MergeArrayElementInitialSections_2 { input: files_to_merge = MergeArrayElementInitialSections_1.merged_file, merged_file_name = "EBR_initial_sections.txt" }
        call Utils.MergeFiles as MergeArrayElementFinalSections_2 { input: files_to_merge = MergeArrayElementFinalSections_1.merged_file, merged_file_name = "EBR_final_sections.txt" }

        # Merge the 10x stats:
        call Utils.MergeCountTsvFiles as Merge10XStats_1 { input: count_tsv_files = AnnotateArrayElements.stats }

        # Merge all Aligned array elements together for this Subread BAM:
        call Utils.MergeBams as MergeAlignedAnnoatatedArrayElementChunk { input: bams = AnnotateArrayElements.output_bam }

        # Merge all Aligned array elements together for this Subread BAM:
        call Utils.MergeBams as MergeAlignedArrayElementChunk { input: bams = AlignArrayElements.aligned_bam }

        # Merge all CCS bams together for this Subread BAM:
        call Utils.MergeBams as MergeChunks { input: bams = CCS.consensus }

        # Merge all CCS bams together for this Subread BAM:
        call Utils.MergeBams as MergeAlignedChunks { input: bams = AlignCCSReads.aligned_bam }

        # Merge all CCS reports together for this Subread BAM:
        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report }

        # Collect metrics on the subreads bam:
        RuntimeAttr subreads_sam_stats_runtime_attrs = object {
            cpu_cores:          2,
            mem_gb:             8,
            disk_gb:            ceil(3 * size(subread_bam, "GiB")),
            boot_disk_gb:       10,
            preemptible_tries:  0,
            max_retries:        1,
            docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.8"
        }
        call AM.SamtoolsStats as CalcSamStatsOnInputBam {
            input:
                bam = subread_bam,
                runtime_attr_override = subreads_sam_stats_runtime_attrs
        }
        call FF.FinalizeToDir as FinalizeSamStatsOnInputBam {
            input:
                # an unfortunate hard-coded path here:
                outdir = outdir + "/" + DIR + "/" + WdlExecutionStartTimestamp.timestamp_string + "/metrics/input_bam_stats",
                files = [
                    CalcSamStatsOnInputBam.raw_stats,
                    CalcSamStatsOnInputBam.summary_stats,
                    CalcSamStatsOnInputBam.first_frag_qual,
                    CalcSamStatsOnInputBam.last_frag_qual,
                    CalcSamStatsOnInputBam.first_frag_gc_content,
                    CalcSamStatsOnInputBam.last_frag_gc_content,
                    CalcSamStatsOnInputBam.acgt_content_per_cycle,
                    CalcSamStatsOnInputBam.insert_size,
                    CalcSamStatsOnInputBam.read_length_dist,
                    CalcSamStatsOnInputBam.indel_distribution,
                    CalcSamStatsOnInputBam.indels_per_cycle,
                    CalcSamStatsOnInputBam.coverage_distribution,
                    CalcSamStatsOnInputBam.gc_depth
                ]
        }
    }

    # Merge all sharded merged sharded outputs of ExtractBoundedReadSectionsTask:
    # Phew.  This is a lot of merging.
    call Utils.MergeFiles as MergeArrayElementSubShards_3 { input: files_to_merge = MergeArrayElementSubShards_2.merged_file, merged_file_name = "EBR_extracted_reads.fasta" }
    call Utils.MergeFiles as MergeArrayElementRejectedReads_3 { input: files_to_merge = MergeArrayElementRejectedReads_2.merged_file, merged_file_name = "EBR_rejected_reads.txt" }
    call Utils.MergeFiles as MergeArrayElementMarkerAlignments_3 { input: files_to_merge = MergeArrayElementMarkerAlignments_2.merged_file, merged_file_name = "EBR_marker_alignments.txt" }
    call Utils.MergeFiles as MergeArrayElementInitialSections_3 { input: files_to_merge = MergeArrayElementInitialSections_2.merged_file, merged_file_name = "EBR_initial_sections.txt" }
    call Utils.MergeFiles as MergeArrayElementFinalSections_3 { input: files_to_merge = MergeArrayElementFinalSections_2.merged_file, merged_file_name = "EBR_final_sections.txt" }

    # Merge the 10x merged stats:
    call Utils.MergeCountTsvFiles as Merge10XStats_2 { input: count_tsv_files = Merge10XStats_1.merged_tsv }

    # Merge all array element bams together for this flowcell:
    call Utils.MergeBams as MergeAnnotatedAlignedArrayElements { input: bams = MergeAlignedAnnoatatedArrayElementChunk.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }

    # Merge all array element bams together for this flowcell:
    call Utils.MergeBams as MergeAlignedArrayElements { input: bams = MergeAlignedArrayElementChunk.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }

    # Merge all aligned CCS bams together for this flowcell:
    call Utils.MergeBams as MergeAllAlignedCCSBams { input: bams = MergeAlignedChunks.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }

    # Merge all CCS bams together for this flowcell:
    call Utils.MergeBams as MergeAllCCSBams { input: bams = MergeChunks.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }

    # Merge all CCS reports together for this flowcell:
    call PB.MergeCCSReports as MergeAllCCSReports { input: reports = MergeCCSReports.report }

    ##########
    # Metrics and plots
    ##########

    String base_out_dir = outdir + "/" + DIR[0] + "/" + WdlExecutionStartTimestamp.timestamp_string
    String metrics_out_dir = base_out_dir + "/metrics"

    # Aligned CCS Metrics:
    call RM.CalculateAndFinalizeReadMetrics as AlignedCCSMetrics {
        input:
            bam_file = MergeAllAlignedCCSBams.merged_bam,
            bam_index = MergeAllAlignedCCSBams.merged_bai,
            ref_dict = ref_fasta_dict,

            base_metrics_out_dir = metrics_out_dir + "/aligned_ccs_metrics"
    }

    # Aligned Array Element Metrics:
    call RM.CalculateAndFinalizeAlternateReadMetrics as AlignedArrayElementMetrics {
        input:
            bam_file = MergeAlignedArrayElements.merged_bam,
            bam_index = MergeAlignedArrayElements.merged_bai,
            ref_dict = ref_fasta_dict,

            base_metrics_out_dir = metrics_out_dir + "/aligned_array_element_metrics"
    }

    # RNASeqMetrics:
    call AM.RnaSeqMetrics as ArrayElementRnaSeqMetrics {
        input:
            bam = MergeAlignedArrayElements.merged_bam,
            bai = MergeAlignedArrayElements.merged_bai,
            ref_flat = ref_flat_annotations
    }

    ##########
    # Create Report:
    ##########

    ## NOTE: This assumes ONE file for both the raw input and the 10x array element stats!
    ##       This should be fixed in version 2.
    call JUPYTER.PB10xMasSeqSingleFlowcellReport as GenerateInteractiveReport {
        input:
            notebook_template              = jupyter_template_interactive,

            subreads_stats                 = CalcSamStatsOnInputBam.raw_stats[0],
            ccs_reads_stats                = AlignedCCSMetrics.sam_stats_raw_stats,
            array_elements_stats           = AlignedArrayElementMetrics.sam_stats_raw_stats,
            ccs_report_file                = MergeAllCCSReports.report,

            ccs_bam_file                   = MergeAllCCSBams.merged_bam,
            array_element_bam_file         = MergeAlignedArrayElements.merged_bam,

            ebr_element_marker_alignments  = MergeArrayElementMarkerAlignments_3.merged_file,
            ebr_initial_section_alignments = MergeArrayElementInitialSections_3.merged_file,
            ebr_final_section_alignments   = MergeArrayElementFinalSections_3.merged_file,
            ebr_bounds_file                = boundaries_file,

            ten_x_metrics_file             = Merge10XStats_2.merged_tsv,
            rna_seq_metrics_file           = ArrayElementRnaSeqMetrics.rna_metrics,

            workflow_dot_file              = workflow_dot_file,
            prefix                         = "interactive"
    }

    ## NOTE: This assumes ONE file for both the raw input and the 10x array element stats!
    ##       This should be fixed in version 2.
    call JUPYTER.PB10xMasSeqSingleFlowcellReport as GenerateStaticReport {
        input:
            notebook_template              = jupyter_template_static,

            subreads_stats                 = CalcSamStatsOnInputBam.raw_stats[0],
            ccs_reads_stats                = AlignedCCSMetrics.sam_stats_raw_stats,
            array_elements_stats           = AlignedArrayElementMetrics.sam_stats_raw_stats,
            ccs_report_file                = MergeAllCCSReports.report,

            ccs_bam_file                   = MergeAllCCSBams.merged_bam,
            array_element_bam_file         = MergeAlignedArrayElements.merged_bam,

            ebr_element_marker_alignments  = MergeArrayElementMarkerAlignments_3.merged_file,
            ebr_initial_section_alignments = MergeArrayElementInitialSections_3.merged_file,
            ebr_final_section_alignments   = MergeArrayElementFinalSections_3.merged_file,
            ebr_bounds_file                = boundaries_file,

            ten_x_metrics_file             = Merge10XStats_2.merged_tsv,
            rna_seq_metrics_file           = ArrayElementRnaSeqMetrics.rna_metrics,

            workflow_dot_file              = workflow_dot_file,
            prefix                         = "static"
    }

    ##########
    # Finalize
    ##########

    # TODO: Should make this iterate through all found bam files, not just the first one.  This seems to be the case for all workflows right now though...

    # Finalize the notebooks:
    String interactive_report_dir = metrics_out_dir + "/report_interactive"
    call FF.FinalizeToDir as FinalizeInteractiveReport {
        input:
            files = [
                GenerateInteractiveReport.populated_notebook,
                GenerateInteractiveReport.html_report,
                GenerateInteractiveReport.pdf_report
            ],
            outdir = interactive_report_dir
    }
    String static_report_dir = metrics_out_dir + "/report_static"
    call FF.FinalizeToDir as FinalizeStaticReport {
        input:
            files = [
                GenerateStaticReport.populated_notebook,
                GenerateStaticReport.html_report,
                GenerateStaticReport.pdf_report
            ],
            outdir = static_report_dir
    }

    # Finalize all the Extracted Bounded Regions data:
    String extractBoundedRegionsDir = base_out_dir + "/extract_bounded_regions"
    call FF.FinalizeToDir as FinalizeEbrData {
        input:
            files = [
                MergeArrayElementSubShards_3.merged_file,
                MergeArrayElementRejectedReads_3.merged_file,
                MergeArrayElementMarkerAlignments_3.merged_file,
                MergeArrayElementInitialSections_3.merged_file,
                MergeArrayElementFinalSections_3.merged_file
            ],
            outdir = extractBoundedRegionsDir
    }

    # Finalize all the 10x metrics here:
    String tenXToolMetricsDir = metrics_out_dir + "/ten_x_tool_metrics"
    scatter ( i in range(length(AnnotateArrayElements.output_bam[0]))) {
        call FF.FinalizeToDir as FinalizeTenXRgStats {
            input:
                files = [
                    AnnotateArrayElements.barcode_stats[0][i],
                    AnnotateArrayElements.starcode[0][i],
                    AnnotateArrayElements.stats[0][i],
                    AnnotateArrayElements.timing_info[0][i]
                ],
                outdir = tenXToolMetricsDir + "/" + i
        }
    }

    call FF.FinalizeToDir as FinalizeAnnotatedArrayElements {
        input:
            files = [ MergeAnnotatedAlignedArrayElements.merged_bam, MergeAnnotatedAlignedArrayElements.merged_bai ],
            outdir = base_out_dir + "/annotated_aligned_array_elements"
    }

    call FF.FinalizeToDir as FinalizeRnaSeqMetrics {
        input:
            files = [ ArrayElementRnaSeqMetrics.rna_metrics ],
            outdir = metrics_out_dir
    }

    call FF.FinalizeToDir as FinalizeArrayElements {
        input:
            files = [ MergeAlignedArrayElements.merged_bam, MergeAlignedArrayElements.merged_bai ],
            outdir = base_out_dir + "/aligned_array_elements"
    }

    call FF.FinalizeToDir as FinalizeAlignedCCSBams {
        input:
            files = [ MergeAllAlignedCCSBams.merged_bam, MergeAllAlignedCCSBams.merged_bai ],
            outdir = base_out_dir + "/merged_bams/aligned"
    }

    call FF.FinalizeToDir as FinalizeCCSBams {
        input:
            files = [ MergeAllCCSBams.merged_bam, MergeAllCCSBams.merged_bai ],
            outdir = base_out_dir + "/merged_bams/unaligned"
    }

    call FF.FinalizeToDir as FinalizeCCSMetrics {
        input:
            files = [ MergeAllCCSReports.report ],
            outdir = metrics_out_dir + "/ccs_metrics"
    }
}
