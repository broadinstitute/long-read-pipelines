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
import "tasks/Longbow.wdl" as LONGBOW

import "tasks/TranscriptAnalysis/UMI_Tools.wdl" as UMI_TOOLS
import "tasks/TranscriptAnalysis/Postprocessing_Tasks.wdl" as TX_POST

workflow MasSeqAlignmentOnly {

    meta {
        description : "This workflow will process and align MAS-seq data.  No UMI / CBC calculations are performed."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String gcs_input_dir
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/MasSeqAlignmentOnly"

        # NOTE: Reference for un-split CCS reads:
        File ref_fasta =  "gs://broad-dsde-methods-long-reads/resources/references/grch38/Homo_sapiens_assembly38.fasta"
        File ref_fasta_index = "gs://broad-dsde-methods-long-reads/resources/references/grch38/Homo_sapiens_assembly38.fasta.fai"
        File ref_fasta_dict = "gs://broad-dsde-methods-long-reads/resources/references/grch38/Homo_sapiens_assembly38.dict"

        # NOTE: Reference for array elements:
        File transcriptome_ref_fasta =  "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.pc_transcripts.fa"
        File transcriptome_ref_fasta_index = "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.pc_transcripts.fa.fai"
        File transcriptome_ref_fasta_dict = "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.pc_transcripts.dict"

        # Default here is 0 because ccs uncorrected reads all seem to have RQ = -1.
        # All pathologically long reads also have RQ = -1.
        # This way we preserve the vast majority of the data, even if it has low quality.
        # We can filter it out at later steps.
        Float min_read_quality = 0.0
        Int max_reclamation_length = 60000

        String mas_seq_model = "mas15"

        String? sample_name
    }

    parameter_meta {
        gcs_input_dir : "Input folder on GCS in which to search for BAM files to process."
        gcs_out_root_dir : "Root output GCS folder in which to place results of this workflow."

        ref_fasta : "FASTA file containing the reference sequence to which the input data should be aligned before splitting into array elements."
        ref_fasta_index : "FASTA index file for the given ref_fasta file."
        ref_fasta_dict : "Sequence dictionary file for the given ref_fasta file."

        transcriptome_ref_fasta : "FASTA file containing the reference sequence to which the array elements should be aligned."
        transcriptome_ref_fasta_index : "FASTA index file for the given transcriptome_ref_fasta file."
        transcriptome_ref_fasta_dict : "Sequence dictionary file for the given transcriptome_ref_fasta file."

        min_read_quality : "[optional] Minimum read quality for reads to have to be included in our data (Default: 0.0)."
        max_reclamation_length : "[optional] Maximum length (in bases) that a read can be to attempt to reclaim from CCS rejection (Default: 60000)."

        sample_name : "[optional] The name of the sample to associate with the data in this workflow."
    }

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_01_WdlExecutionStartTimestamp { input: }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams as t_02_FindBams { input: gcs_input_dir = gcs_input_dir }

    # Check here if we found ccs bams or subread bams:
    Boolean use_subreads = t_02_FindBams.has_subreads

    # Make sure we have **EXACTLY** one bam file to run on:
    if (length(t_02_FindBams.ccs_bams) != 1) {
        call Utils.FailWithWarning as t_03_FailOnMultiBamFiles { input: warning = "Error: Multiple BAM files found.  Cannot continue!" }
     }

    # Alias our bam file so we can work with it easier:
    File reads_bam = t_02_FindBams.ccs_bams[0]

    call PB.GetRunInfo as t_04_GetRunInfo { input: subread_bam = reads_bam }

    String SM  = select_first([sample_name, t_04_GetRunInfo.run_info["SM"]])
    String PL  = "PACBIO"
    String PU  = t_04_GetRunInfo.run_info["PU"]
    String DT  = t_04_GetRunInfo.run_info["DT"]
    String ID  = PU
    String DS  = t_04_GetRunInfo.run_info["DS"]
    String DIR = SM + "." + ID

    File read_pbi = sub(reads_bam, ".bam$", ".bam.pbi")
    call PB.ShardLongReads as t_05_ShardLongReads {
        input:
            unaligned_bam = reads_bam,
            unaligned_pbi = read_pbi,
            prefix = SM + "_shard",
            num_shards = 50,
    }

    scatter (sharded_reads in t_05_ShardLongReads.unmapped_shards) {

        ## No more preemption on this sharding - takes too long otherwise.
        RuntimeAttr disable_preemption_runtime_attrs = object {
            preemptible_tries: 0
        }

        String fbmrq_prefix = basename(sharded_reads, ".bam")

        # Filter out the kinetics tags from PB files:
        call PB.RemoveKineticsTags as t_06_RemoveKineticsTags {
            input:
                bam = sharded_reads,
                prefix = SM + "_kinetics_removed"
        }

        # Handle setting up the things that we need for further processing of CCS-only reads:
        call PB.FindCCSReport as t_07_FindCCSReport {
            input:
                gcs_input_dir = gcs_input_dir
        }

        # 1 - filter the reads by the minimum read quality:
        call Utils.Bamtools as t_08_FilterByMinQual {
            input:
                bamfile = t_06_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_good_reads",
                cmd = "filter",
                args = '-tag "rq":">=' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        call LONGBOW.Annotate as t_09_LongbowAnnotateReads {
            input:
                reads = t_08_FilterByMinQual.bam_out,
                model = mas_seq_model
        }

        call PB.PBIndex as t_10_PbIndexLongbowAnnotatedReads {
            input:
                bam = t_09_LongbowAnnotateReads.annotated_bam
        }

        # Shard these reads even wider so we can make sure we don't run out of memory:
        call PB.ShardLongReads as t_11_ShardCorrectedReads {
            input:
                unaligned_bam = t_09_LongbowAnnotateReads.annotated_bam,
                unaligned_pbi = t_10_PbIndexLongbowAnnotatedReads.pbindex,
                prefix = SM + "_longbow_annotated_subshard",
                num_shards = 10,
        }

        # Segment our arrays into individual array elements:
        scatter (corrected_shard in t_11_ShardCorrectedReads.unmapped_shards) {
            call LONGBOW.Segment as t_12_SegmentAnnotatedReads {
                input:
                    annotated_reads = corrected_shard,
                    model = mas_seq_model
            }
        }

        # Merge all outputs of Longbow Annotate / Segment:
        call Utils.MergeBams as t_13_MergeArrayElements_1 {
            input:
                bams = t_12_SegmentAnnotatedReads.segmented_bam,
                prefix = SM + "_ArrayElements_intermediate_1"
        }

        # Align our array elements:
        call AR.Minimap2 as t_14_AlignArrayElements {
            input:
                reads      = [ t_13_MergeArrayElements_1.merged_bam ],
                ref_fasta  = transcriptome_ref_fasta,
                map_preset = "splice:hq"
        }

        call AR.Minimap2 as t_15_AlignArrayElementsToGenome {
            input:
                reads      = [ t_13_MergeArrayElements_1.merged_bam ],
                ref_fasta  = ref_fasta,
                map_preset = "splice:hq"
        }

        # We need to restore the annotations we created with the 10x tool to the aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_16_RestoreAnnotationsToTranscriptomeAlignedBam {
            input:
                annotated_bam_file = t_13_MergeArrayElements_1.merged_bam,
                aligned_bam_file = t_14_AlignArrayElements.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # We need to restore the annotations we created with the 10x tool to the aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_17_RestoreAnnotationsToGenomeAlignedBam {
            input:
                annotated_bam_file = t_13_MergeArrayElements_1.merged_bam,
                aligned_bam_file = t_15_AlignArrayElementsToGenome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # To properly count our transcripts we must throw away the non-primary and unaligned reads:
        RuntimeAttr filterReadsAttrs = object {
            cpu_cores: 4,
            preemptible_tries: 0
        }
        call Utils.FilterReadsBySamFlags as t_18_RemoveUnmappedAndNonPrimaryTranscriptomeReads {
            input:
                bam = t_16_RestoreAnnotationsToTranscriptomeAlignedBam.output_bam,
                sam_flags = "2308",
                prefix = SM + "_ArrayElements_Annotated_Transcriptome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }
        call Utils.FilterReadsBySamFlags as t_19_RemoveUnmappedAndNonPrimaryGenomeReads {
            input:
                bam = t_17_RestoreAnnotationsToGenomeAlignedBam.output_bam,
                sam_flags = "2308",
                prefix = SM + "_ArrayElements_Annotated_Genome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }
    }

    #################################################
    #   __  __
    #  |  \/  | ___ _ __ __ _  ___
    #  | |\/| |/ _ \ '__/ _` |/ _ \
    #  | |  | |  __/ | | (_| |  __/
    #  |_|  |_|\___|_|  \__, |\___|
    #                   |___/
    #
    # Merge all the sharded files we created above into files for this
    # input bam file.
    #################################################

    # Sequel IIe Data.
    # CCS Passed:
    call Utils.MergeBams as t_20_MergeCCSRqFilteredReads { input: bams = t_08_FilterByMinQual.bam_out, prefix = SM + "_ccs_reads" }
    call Utils.MergeBams as t_21_MergeAnnotatedCCSReads { input: bams = t_09_LongbowAnnotateReads.annotated_bam, prefix = SM + "_ccs_reads_annotated" }

    # Merge all CCS bams together for this Subread BAM:
    RuntimeAttr merge_extra_cpu_attrs = object {
        cpu_cores: 4
    }
    call Utils.MergeBams as t_22_MergeRawArrayElements { input: bams = t_13_MergeArrayElements_1.merged_bam, prefix = SM + "_raw_array_elements" }
    call Utils.MergeBams as t_23_MergeTranscriptomeAlignedArrayElements { input: bams = t_16_RestoreAnnotationsToTranscriptomeAlignedBam.output_bam, prefix = SM + "_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_24_MergeGenomeAlignedArrayElements { input: bams = t_17_RestoreAnnotationsToGenomeAlignedBam.output_bam, prefix = SM + "_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_25_MergePrimaryTranscriptomeAlignedArrayElements { input: bams = t_18_RemoveUnmappedAndNonPrimaryTranscriptomeReads.output_bam, prefix = SM + "_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_26_MergePrimaryGenomeAlignedArrayElements { input: bams = t_19_RemoveUnmappedAndNonPrimaryGenomeReads.output_bam, prefix = SM + "_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    #################################################
    #   ___      ____
    #  / _ \    / ___|
    # | | | |  | |
    # | |_| |  | |___
    #  \__\_\   \____|
    #
    #################################################

    call LONGBOW.Stats as t_27_longbow_stats {
        input:
            reads = t_21_MergeAnnotatedCCSReads.merged_bam,
            model = mas_seq_model,
            prefix = SM
    }

    ######################################################################
    #             _____ _             _ _
    #            |  ___(_)_ __   __ _| (_)_______
    #            | |_  | | '_ \ / _` | | |_  / _ \
    #            |  _| | | | | | (_| | | |/ /  __/
    #            |_|   |_|_| |_|\__,_|_|_/___\___|
    #
    ######################################################################

    # NOTE: We key all finalization steps on the static report.
    #       This will prevent incomplete runs from being placed in the output folders.
    String base_out_dir = outdir + "/" + DIR + "/" + t_01_WdlExecutionStartTimestamp.timestamp_string

    File final_ccs_report = t_07_FindCCSReport.ccs_report[0]

    String array_element_dir = base_out_dir + "/array_elements"
    String arrays_dir = base_out_dir + "/array_reads"

    ##############################################################################################################
    # Finalize the final annotated, aligned array elements:
    call FF.FinalizeToDir as t_28_FinalizeCCSReads {
        input:
            files = [
                t_20_MergeCCSRqFilteredReads.merged_bam,
                t_20_MergeCCSRqFilteredReads.merged_bai,
                t_21_MergeAnnotatedCCSReads.merged_bam,
                t_21_MergeAnnotatedCCSReads.merged_bai,
            ],
            outdir = arrays_dir,
            keyfile = t_27_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the intermediate reads files (from raw CCS corrected reads through split array elements)
    call FF.FinalizeToDir as t_29_FinalizeArrayElements {
        input:
            files = [
                t_22_MergeRawArrayElements.merged_bam,
                t_22_MergeRawArrayElements.merged_bai,
                t_23_MergeTranscriptomeAlignedArrayElements.merged_bam,
                t_23_MergeTranscriptomeAlignedArrayElements.merged_bai,
                t_24_MergeGenomeAlignedArrayElements.merged_bam,
                t_24_MergeGenomeAlignedArrayElements.merged_bai,
                t_25_MergePrimaryTranscriptomeAlignedArrayElements.merged_bam,
                t_25_MergePrimaryTranscriptomeAlignedArrayElements.merged_bai,
                t_26_MergePrimaryGenomeAlignedArrayElements.merged_bam,
                t_26_MergePrimaryGenomeAlignedArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = t_27_longbow_stats.summary_stats
    }

    call FF.FinalizeToDir as t_30_FinalizeCCSReport {
        input:
            files = [
                final_ccs_report
            ],
            outdir = base_out_dir + "/",
            keyfile = t_27_longbow_stats.summary_stats
    }

    call FF.FinalizeToDir as t_31_FinalizeLongbowStats {
        input:
            files = [
                t_27_longbow_stats.summary_stats,
                t_27_longbow_stats.array_length_counts_plot_png,
                t_27_longbow_stats.array_length_counts_plot_svg,
                t_27_longbow_stats.ligation_heatmap_nn_png,
                t_27_longbow_stats.ligation_heatmap_nn_svg,
                t_27_longbow_stats.ligation_heatmap_png,
                t_27_longbow_stats.ligation_heatmap_svg,
                t_27_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_27_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_27_longbow_stats.ligation_heatmap_reduced_png,
                t_27_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/",
            keyfile = t_27_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good:
    call FF.WriteCompletionFile as t_32_WriteCompletionFile {
        input:
            outdir = base_out_dir + "/",
            keyfile =  t_27_longbow_stats.summary_stats
    }
}
