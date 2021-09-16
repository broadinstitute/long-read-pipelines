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

        # Set up some meta parameters here so we can adjust for when we want things to go VERY fast:
        Int primary_scatter_width = 50
        Int secondary_scatter_width = 10
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

        primary_scatter_width : "[optional] Width to use for the primary scatter operation on this dataset (default: 50)."
        secondary_scatter_width : "[optional] Width to use for the secondary (nested) scatter operation on this dataset (default: 10)."
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

    # Define some attributes for later:
    RuntimeAttr disable_preemption_runtime_attrs = object {
        preemptible_tries: 0
    }
    RuntimeAttr filterReadsAttrs = object {
        cpu_cores: 4,
        preemptible_tries: 0
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
            num_shards = primary_scatter_width,
    }

    call PB.FindCCSReport as t_06_FindCCSReport {
        input:
            gcs_input_dir = gcs_input_dir
    }

    scatter (sharded_reads in t_05_ShardLongReads.unmapped_shards) {

        #######################################################################
        ## Setup
        #######################################################################

        String fbmrq_prefix = basename(sharded_reads, ".bam")

        # Filter out the kinetics tags from PB files:
        call PB.RemoveKineticsTags as t_07_RemoveKineticsTags {
            input:
                bam = sharded_reads,
                prefix = SM + "_kinetics_removed"
        }

        #######################################################################
        ## Handle CCS Corrected Reads
        #######################################################################

        # 1 - filter the reads by the minimum read quality:
        call Utils.Bamtools as t_08_FilterByMinQual {
            input:
                bamfile = t_07_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_good_reads",
                cmd = "filter",
                args = '-tag "rq":">=' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # Annotate our CCS Corrected reads:
        call LONGBOW.Annotate as t_09_LongbowAnnotateCCSReads {
            input:
                reads = t_08_FilterByMinQual.bam_out,
                model = mas_seq_model
        }

        # 6: Longbow filter ccs reclaimable reads
        call LONGBOW.Filter as t_10_FilterCCSReads {
            input:
                bam = t_09_LongbowAnnotateCCSReads.annotated_bam,
                prefix = SM + "_ccs_corrected_subshard",
                model = mas_seq_model
        }

        call PB.PBIndex as t_11_PbIndexLongbowAnnotatedCCSPassedReads {
            input:
                bam = t_10_FilterCCSReads.passed_reads
        }

        # Shard these reads even wider so we can make sure we don't run out of memory:
        call PB.ShardLongReads as t_12_ShardCorrectedReads {
            input:
                unaligned_bam = t_10_FilterCCSReads.passed_reads,
                unaligned_pbi = t_11_PbIndexLongbowAnnotatedCCSPassedReads.pbindex,
                prefix = SM + "_ccs_corrected_longbow_annotated_subshard",
                num_shards = secondary_scatter_width,
        }

        # Segment our arrays into individual array elements:
        scatter (corrected_shard in t_12_ShardCorrectedReads.unmapped_shards) {
            call LONGBOW.Segment as t_13_SegmentCCSAnnotatedReads {
                input:
                    annotated_reads = corrected_shard,
                    model = mas_seq_model
            }
        }

        # Merge all outputs of Longbow Annotate / Segment:
        call Utils.MergeBams as t_14_MergeCCSArrayElements_1 {
            input:
                bams = t_13_SegmentCCSAnnotatedReads.segmented_bam,
                prefix = SM + "_ccs_corrected_ArrayElements_intermediate_1"
        }

        # Align our array elements to the transcriptome:
        call AR.Minimap2 as t_15_AlignCCSArrayElements {
            input:
                reads      = [ t_14_MergeCCSArrayElements_1.merged_bam ],
                ref_fasta  = transcriptome_ref_fasta,
                map_preset = "asm20"
        }

        # Align our array elements to the genome in splice-aware mode:
        call AR.Minimap2 as t_16_AlignCCSArrayElementsToGenome {
            input:
                reads      = [ t_14_MergeCCSArrayElements_1.merged_bam ],
                ref_fasta  = ref_fasta,
                map_preset = "splice:hq"
        }

        # We need to restore the annotations we created with the 10x tool to the transcriptome aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_17_RestoreAnnotationsToTranscriptomeAlignedCCSBam {
            input:
                annotated_bam_file = t_14_MergeCCSArrayElements_1.merged_bam,
                aligned_bam_file = t_15_AlignCCSArrayElements.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # We need to restore the annotations we created with the 10x tool to the genome aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_18_RestoreAnnotationsToGenomeAlignedCCSBam {
            input:
                annotated_bam_file = t_14_MergeCCSArrayElements_1.merged_bam,
                aligned_bam_file = t_16_AlignCCSArrayElementsToGenome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # To properly count our transcripts we must throw away the non-primary and unaligned reads:
        # Filter out all non-primary transcriptome alignments:
        call Utils.FilterReadsBySamFlags as t_19_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads {
            input:
                bam = t_17_RestoreAnnotationsToTranscriptomeAlignedCCSBam.output_bam,
                sam_flags = "2308",
                prefix = SM + "_ccs_corrected_ArrayElements_Annotated_Transcriptome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }
        # Filter out all non-primary genome alignments:
        call Utils.FilterReadsBySamFlags as t_20_RemoveUnmappedAndNonPrimaryCCSGenomeReads {
            input:
                bam = t_18_RestoreAnnotationsToGenomeAlignedCCSBam.output_bam,
                sam_flags = "2308",
                prefix = SM + "_ccs_corrected_ArrayElements_Annotated_Genome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }

        #######################################################################
        ## Handle CCS Uncorrected Reads
        #######################################################################

        # 1.5 - Get the "rejected" reads:
        call Utils.Bamtools as t_21_GetCcsRejectedReads {
            input:
                bamfile = t_07_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_rejected_reads",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # 2 - Get reads we can reclaim:
        call Utils.Bamtools as t_22_ExtractCcsReclaimableReads {
            input:
                bamfile = t_07_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_reads_for_ccs_reclamation",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '" -length "<=' + max_reclamation_length + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # Annotate our CCS uncorrected (reclaimable) reads
        call LONGBOW.Annotate as t_23_AnnotateReclaimableReads {
            input:
                reads = t_22_ExtractCcsReclaimableReads.bam_out,
                model = mas_seq_model
        }

        call PB.PBIndex as t_24_PbIndexLongbowAnnotatedReclaimedReads {
            input:
                bam = t_23_AnnotateReclaimableReads.annotated_bam
        }

        # 6: Longbow filter ccs reclaimable reads
        call LONGBOW.Filter as t_25_FilterReclaimableReads {
            input:
                bam = t_23_AnnotateReclaimableReads.annotated_bam,
                bam_pbi = t_24_PbIndexLongbowAnnotatedReclaimedReads.pbindex,
                prefix = SM + "_subshard",
                model = mas_seq_model
        }

        call PB.PBIndex as t_26_PbIndexLongbowAnnotatedReclaimedPassedReads {
            input:
                bam = t_25_FilterReclaimableReads.passed_reads
        }

        # Shard these reads even wider so we can make sure we don't run out of memory:
        call PB.ShardLongReads as t_27_ShardReclaimedReads {
            input:
                unaligned_bam = t_25_FilterReclaimableReads.passed_reads,
                unaligned_pbi = t_26_PbIndexLongbowAnnotatedReclaimedPassedReads.pbindex,
                prefix = SM + "_ccs_reclaimed_longbow_annotated_subshard",
                num_shards = secondary_scatter_width,
        }

        # Segment our arrays into individual array elements:
        scatter (corrected_shard in t_27_ShardReclaimedReads.unmapped_shards) {
            call LONGBOW.Segment as t_28_SegmentReclaimedAnnotatedReads {
                input:
                    annotated_reads = corrected_shard,
                    model = mas_seq_model
            }
        }

        # Merge all outputs of Longbow Annotate / Segment:
        call Utils.MergeBams as t_29_MergeReclaimedArrayElements_1 {
            input:
                bams = t_28_SegmentReclaimedAnnotatedReads.segmented_bam,
                prefix = SM + "_ccs_reclaimed_ArrayElements_intermediate_1"
        }

        # Align our array elements to the transcriptome:
        call AR.Minimap2 as t_30_AlignReclaimedArrayElements {
            input:
                reads      = [ t_29_MergeReclaimedArrayElements_1.merged_bam ],
                ref_fasta  = transcriptome_ref_fasta,
                map_preset = "asm20"
        }

        # Align our array elements to the genome in splice-aware mode:
        call AR.Minimap2 as t_31_AlignReclaimedArrayElementsToGenome {
            input:
                reads      = [ t_29_MergeReclaimedArrayElements_1.merged_bam ],
                ref_fasta  = ref_fasta,
                map_preset = "splice:hq"
        }

        # We need to restore the annotations we created with the 10x tool to the transcriptome aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_32_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam {
            input:
                annotated_bam_file = t_29_MergeReclaimedArrayElements_1.merged_bam,
                aligned_bam_file = t_30_AlignReclaimedArrayElements.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # We need to restore the annotations we created with the 10x tool to the genome aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_33_RestoreAnnotationsToGenomeAlignedReclaimedBam {
            input:
                annotated_bam_file = t_29_MergeReclaimedArrayElements_1.merged_bam,
                aligned_bam_file = t_31_AlignReclaimedArrayElementsToGenome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # To properly count our transcripts we must throw away the non-primary and unaligned reads:
        # Filter out all non-primary transcriptome alignments:
        call Utils.FilterReadsBySamFlags as t_34_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads {
            input:
                bam = t_32_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam.output_bam,
                sam_flags = "2308",
                prefix = SM + "_ccs_reclaimed_ArrayElements_Annotated_Transcriptome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }
        # Filter out all non-primary genome alignments:
        call Utils.FilterReadsBySamFlags as t_35_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads {
            input:
                bam = t_33_RestoreAnnotationsToGenomeAlignedReclaimedBam.output_bam,
                sam_flags = "2308",
                prefix = SM + "_ccs_reclaimed_ArrayElements_Annotated_Genome_Aligned_PrimaryOnly",
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

    RuntimeAttr merge_extra_cpu_attrs = object {
        cpu_cores: 4
    }

    #################################################
    # CCS Passed:
    call Utils.MergeBams as t_36_MergeCCSReads { input: bams = t_08_FilterByMinQual.bam_out, prefix = SM + "_ccs_reads" }
    call Utils.MergeBams as t_37_MergeAnnotatedCCSReads { input: bams = t_09_LongbowAnnotateCCSReads.annotated_bam, prefix = SM + "_ccs_reads_annotated" }

    # Merge all CCS bams together for this Subread BAM:
    call Utils.MergeBams as t_38_MergePassedCCSReads { input: bams = t_10_FilterCCSReads.passed_reads, prefix = SM + "_ccs_longbow_passed_reads" }
    call Utils.MergeBams as t_39_MergeFailedCCSReads { input: bams = t_10_FilterCCSReads.failed_reads, prefix = SM + "_ccs_longbow_failed_reads" }
    call Utils.MergeBams as t_40_MergeCCSArrayElements { input: bams = t_14_MergeCCSArrayElements_1.merged_bam, prefix = SM + "_ccs_array_elements" }
    call Utils.MergeBams as t_41_MergeTranscriptomeAlignedCCSArrayElements { input: bams = t_17_RestoreAnnotationsToTranscriptomeAlignedCCSBam.output_bam, prefix = SM + "_ccs_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_42_MergeGenomeAlignedCCSArrayElements { input: bams = t_18_RestoreAnnotationsToGenomeAlignedCCSBam.output_bam, prefix = SM + "_ccs_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_43_MergePrimaryTranscriptomeAlignedCCSArrayElements { input: bams = t_19_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads.output_bam, prefix = SM + "_ccs_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_44_MergePrimaryGenomeAlignedCCSArrayElements { input: bams = t_20_RemoveUnmappedAndNonPrimaryCCSGenomeReads.output_bam, prefix = SM + "_ccs_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    #################################################
    # CCS Reclaimed:
    call Utils.MergeBams as t_45_MergeCCSReclaimableReads { input: bams = t_22_ExtractCcsReclaimableReads.bam_out, prefix = SM + "_reclaimable_reads" }
    call Utils.MergeBams as t_46_MergeAnnotatedCCSReclaimableReads { input: bams = t_23_AnnotateReclaimableReads.annotated_bam, prefix = SM + "_reclaimable_reads_annotated" }

    # Merge all CCS bams together for this Subread BAM:
    call Utils.MergeBams as t_47_MergePassedCCSReclaimedReads { input: bams = t_25_FilterReclaimableReads.passed_reads, prefix = SM + "_reclaimable_longbow_passed_reads" }
    call Utils.MergeBams as t_48_MergeFailedCCSReclaimableReads { input: bams = t_25_FilterReclaimableReads.failed_reads, prefix = SM + "_reclaimed_longbow_failed_reads" }
    call Utils.MergeBams as t_49_MergeCCSReclaimedArrayElements { input: bams = t_29_MergeReclaimedArrayElements_1.merged_bam, prefix = SM + "_reclaimed_array_elements" }
    call Utils.MergeBams as t_50_MergeTranscriptomeAlignedCCSReclaimedArrayElements { input: bams = t_32_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam.output_bam, prefix = SM + "_reclaimed_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_51_MergeGenomeAlignedCCSReclaimedArrayElements { input: bams = t_33_RestoreAnnotationsToGenomeAlignedReclaimedBam.output_bam, prefix = SM + "_reclaimed_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_52_MergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements { input: bams = t_34_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads.output_bam, prefix = SM + "_reclaimed_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_53_MergePrimaryGenomeAlignedCCSReclaimedArrayElements { input: bams = t_35_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads.output_bam, prefix = SM + "_reclaimed_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    #################################################
    # Overall:
    call Utils.MergeBams as t_54_MergeAllAnnotatedReads { input: bams = flatten([t_09_LongbowAnnotateCCSReads.annotated_bam, t_23_AnnotateReclaimableReads.annotated_bam]), prefix = SM + "_all_annotated_reads" }
    call Utils.MergeBams as t_55_MergeAllLongbowPassedReads { input: bams = flatten([t_10_FilterCCSReads.passed_reads, t_25_FilterReclaimableReads.passed_reads]), prefix = SM + "_all_longbow_passed_reads" }
    call Utils.MergeBams as t_56_MergeAllLongbowFailedReads { input: bams = flatten([t_10_FilterCCSReads.failed_reads, t_25_FilterReclaimableReads.failed_reads]), prefix = SM + "_all_longbow_failed_reads" }

    call Utils.MergeBams as t_57_MergeAllArrayElements { input: bams = flatten([t_14_MergeCCSArrayElements_1.merged_bam, t_29_MergeReclaimedArrayElements_1.merged_bam]), prefix = SM + "_reclaimed_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_58_MergeAllTranscriptomeAlignedArrayElements { input: bams = t_32_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam.output_bam, prefix = SM + "_reclaimed_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_59_MergeAllGenomeAlignedArrayElements { input: bams = t_33_RestoreAnnotationsToGenomeAlignedReclaimedBam.output_bam, prefix = SM + "_reclaimed_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_60_MergeAllPrimaryTranscriptomeAlignedArrayElements { input: bams = t_34_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads.output_bam, prefix = SM + "_reclaimed_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_61_MergeAllPrimaryGenomeAlignedArrayElements { input: bams = t_35_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads.output_bam, prefix = SM + "_reclaimed_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }


    #################################################
    #   ___      ____
    #  / _ \    / ___|
    # | | | |  | |
    # | |_| |  | |___
    #  \__\_\   \____|
    #
    #################################################

    # Get stats on CCS reads:
    call LONGBOW.Stats as t_62_CCS_longbow_stats {
        input:
            reads = t_37_MergeAnnotatedCCSReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_CCS_Corrected"
    }

    # Get stats on Reclaimable reads:
    call LONGBOW.Stats as t_63_Reclaimable_longbow_stats {
        input:
            reads = t_46_MergeAnnotatedCCSReclaimableReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_CCS_Reclaimable"
    }

    # Get stats on Reclaimed reads:
    call LONGBOW.Stats as t_64_Reclaimed_longbow_stats {
        input:
            reads = t_47_MergePassedCCSReclaimedReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_CCS_Reclaimed"
    }

    # Get stats on All reads (overall stats):
    call LONGBOW.Stats as t_65_Passed_longbow_stats {
        input:
            reads = t_55_MergeAllLongbowPassedReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_All_Longbow_Passed"
    }

    # Get stats on All reads (overall stats):
    call LONGBOW.Stats as t_66_Overall_longbow_stats {
        input:
            reads = t_54_MergeAllAnnotatedReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_Overall"
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

    File final_ccs_report = t_06_FindCCSReport.ccs_report

    String array_element_dir = base_out_dir + "/array_elements"
    String arrays_dir = base_out_dir + "/array_reads"

    ##############################################################################################################
    # Finalize the CCS Array Reads
    call FF.FinalizeToDir as t_67_FinalizeCCSReads {
        input:
            files = [
                t_36_MergeCCSReads.merged_bam,
                t_36_MergeCCSReads.merged_bai,
                t_37_MergeAnnotatedCCSReads.merged_bam,
                t_37_MergeAnnotatedCCSReads.merged_bai,
            ],
            outdir = arrays_dir,
            keyfile = t_66_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Reclaimed Array Reads
    call FF.FinalizeToDir as t_68_FinalizeCCSReclaimedReads {
        input:
            files = [
                t_45_MergeCCSReclaimableReads.merged_bam,
                t_45_MergeCCSReclaimableReads.merged_bai,
                t_46_MergeAnnotatedCCSReclaimableReads.merged_bam,
                t_46_MergeAnnotatedCCSReclaimableReads.merged_bai,
                t_47_MergePassedCCSReclaimedReads.merged_bam,
                t_47_MergePassedCCSReclaimedReads.merged_bai,
                t_48_MergeFailedCCSReclaimableReads.merged_bam,
                t_48_MergeFailedCCSReclaimableReads.merged_bai,
            ],
            outdir = arrays_dir,
            keyfile = t_66_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Overall Array Reads
    call FF.FinalizeToDir as t_69_FinalizeOverallCombinedReads {
        input:
            files = [
                t_54_MergeAllAnnotatedReads.merged_bam,
                t_54_MergeAllAnnotatedReads.merged_bai,
                t_55_MergeAllLongbowPassedReads.merged_bam,
                t_55_MergeAllLongbowPassedReads.merged_bai,
                t_56_MergeAllLongbowFailedReads.merged_bam,
                t_56_MergeAllLongbowFailedReads.merged_bai,
            ],
            outdir = arrays_dir,
            keyfile = t_66_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the CCS Array Element files
    call FF.FinalizeToDir as t_70_FinalizeCCSArrayElements {
        input:
            files = [
                t_40_MergeCCSArrayElements.merged_bam,
                t_40_MergeCCSArrayElements.merged_bai,
                t_41_MergeTranscriptomeAlignedCCSArrayElements.merged_bam,
                t_41_MergeTranscriptomeAlignedCCSArrayElements.merged_bai,
                t_42_MergeGenomeAlignedCCSArrayElements.merged_bam,
                t_42_MergeGenomeAlignedCCSArrayElements.merged_bai,
                t_43_MergePrimaryTranscriptomeAlignedCCSArrayElements.merged_bam,
                t_43_MergePrimaryTranscriptomeAlignedCCSArrayElements.merged_bai,
                t_44_MergePrimaryGenomeAlignedCCSArrayElements.merged_bam,
                t_44_MergePrimaryGenomeAlignedCCSArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = t_66_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Reclaimed Array Element files
    call FF.FinalizeToDir as t_71_FinalizeCCSRecaimedArrayElements {
        input:
            files = [
                t_49_MergeCCSReclaimedArrayElements.merged_bam,
                t_49_MergeCCSReclaimedArrayElements.merged_bai,
                t_50_MergeTranscriptomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_50_MergeTranscriptomeAlignedCCSReclaimedArrayElements.merged_bai,
                t_51_MergeGenomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_51_MergeGenomeAlignedCCSReclaimedArrayElements.merged_bai,
                t_52_MergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_52_MergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements.merged_bai,
                t_53_MergePrimaryGenomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_53_MergePrimaryGenomeAlignedCCSReclaimedArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = t_66_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Overall Array Element files
    call FF.FinalizeToDir as t_72_FinalizeOverallArrayElements {
        input:
            files = [
                t_57_MergeAllArrayElements.merged_bam,
                t_57_MergeAllArrayElements.merged_bai,
                t_58_MergeAllTranscriptomeAlignedArrayElements.merged_bam,
                t_58_MergeAllTranscriptomeAlignedArrayElements.merged_bai,
                t_59_MergeAllGenomeAlignedArrayElements.merged_bam,
                t_59_MergeAllGenomeAlignedArrayElements.merged_bai,
                t_60_MergeAllPrimaryTranscriptomeAlignedArrayElements.merged_bam,
                t_60_MergeAllPrimaryTranscriptomeAlignedArrayElements.merged_bai,
                t_61_MergeAllPrimaryGenomeAlignedArrayElements.merged_bam,
                t_61_MergeAllPrimaryGenomeAlignedArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = t_66_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the CCS report
    call FF.FinalizeToDir as t_73_FinalizeCCSReport {
        input:
            files = [
                final_ccs_report
            ],
            outdir = base_out_dir + "/",
            keyfile = t_66_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the different stats files to separate directories
    call FF.FinalizeToDir as t_74_FinalizeCCSLongbowStats {
        input:
            files = [
                t_62_CCS_longbow_stats.summary_stats,
                t_62_CCS_longbow_stats.array_length_counts_plot_png,
                t_62_CCS_longbow_stats.array_length_counts_plot_svg,
                t_62_CCS_longbow_stats.ligation_heatmap_nn_png,
                t_62_CCS_longbow_stats.ligation_heatmap_nn_svg,
                t_62_CCS_longbow_stats.ligation_heatmap_png,
                t_62_CCS_longbow_stats.ligation_heatmap_svg,
                t_62_CCS_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_62_CCS_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_62_CCS_longbow_stats.ligation_heatmap_reduced_png,
                t_62_CCS_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/CCS_Corrected/",
            keyfile = t_66_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_75_FinalizeReclaimableLongbowStats {
        input:
            files = [
                t_63_Reclaimable_longbow_stats.summary_stats,
                t_63_Reclaimable_longbow_stats.array_length_counts_plot_png,
                t_63_Reclaimable_longbow_stats.array_length_counts_plot_svg,
                t_63_Reclaimable_longbow_stats.ligation_heatmap_nn_png,
                t_63_Reclaimable_longbow_stats.ligation_heatmap_nn_svg,
                t_63_Reclaimable_longbow_stats.ligation_heatmap_png,
                t_63_Reclaimable_longbow_stats.ligation_heatmap_svg,
                t_63_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_63_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_63_Reclaimable_longbow_stats.ligation_heatmap_reduced_png,
                t_63_Reclaimable_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/CCS_Reclaimable/",
            keyfile = t_66_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_76_FinalizeReclaimedLongbowStats {
        input:
            files = [
                t_64_Reclaimed_longbow_stats.summary_stats,
                t_64_Reclaimed_longbow_stats.array_length_counts_plot_png,
                t_64_Reclaimed_longbow_stats.array_length_counts_plot_svg,
                t_64_Reclaimed_longbow_stats.ligation_heatmap_nn_png,
                t_64_Reclaimed_longbow_stats.ligation_heatmap_nn_svg,
                t_64_Reclaimed_longbow_stats.ligation_heatmap_png,
                t_64_Reclaimed_longbow_stats.ligation_heatmap_svg,
                t_64_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_64_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_64_Reclaimed_longbow_stats.ligation_heatmap_reduced_png,
                t_64_Reclaimed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/CCS_Reclaimed/",
            keyfile = t_66_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_77_FinalizeOverallLongbowStats {
        input:
            files = [
                t_66_Overall_longbow_stats.summary_stats,
                t_66_Overall_longbow_stats.array_length_counts_plot_png,
                t_66_Overall_longbow_stats.array_length_counts_plot_svg,
                t_66_Overall_longbow_stats.ligation_heatmap_nn_png,
                t_66_Overall_longbow_stats.ligation_heatmap_nn_svg,
                t_66_Overall_longbow_stats.ligation_heatmap_png,
                t_66_Overall_longbow_stats.ligation_heatmap_svg,
                t_66_Overall_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_66_Overall_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_66_Overall_longbow_stats.ligation_heatmap_reduced_png,
                t_66_Overall_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/Overall/",
            keyfile = t_66_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_78_FinalizeAllPassedLongbowStats {
        input:
            files = [
                t_65_Passed_longbow_stats.summary_stats,
                t_65_Passed_longbow_stats.array_length_counts_plot_png,
                t_65_Passed_longbow_stats.array_length_counts_plot_svg,
                t_65_Passed_longbow_stats.ligation_heatmap_nn_png,
                t_65_Passed_longbow_stats.ligation_heatmap_nn_svg,
                t_65_Passed_longbow_stats.ligation_heatmap_png,
                t_65_Passed_longbow_stats.ligation_heatmap_svg,
                t_65_Passed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_65_Passed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_65_Passed_longbow_stats.ligation_heatmap_reduced_png,
                t_65_Passed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/All_Longbow_Passed/",
            keyfile = t_66_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good
    call FF.WriteCompletionFile as t_79_WriteCompletionFile {
        input:
            outdir = base_out_dir + "/",
            keyfile =  t_66_Overall_longbow_stats.summary_stats
    }
}
