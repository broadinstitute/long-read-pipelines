version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF
import "tasks/AlignReads.wdl" as AR
import "tasks/Ten_X_Tool.wdl" as TENX
import "tasks/Longbow.wdl" as LONGBOW

workflow MasSeqAlignmentOnlyModPrepTests {

    meta {
        description : "This workflow will process and align MAS-seq data.  No UMI / CBC calculations are performed."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File bam_file
        String sample_name
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/MasSeqAlignmentOnlyModPrepTests"

        # NOTE: Reference for un-split CCS reads:
        File ref_fasta =  "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
        File ref_fasta_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
        File ref_fasta_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"

        # NOTE: Reference for array elements:
        File transcriptome_ref_fasta =  "gs://broad-dsde-methods-long-reads-public/resources/gencode_v37/gencode.v37.pc_transcripts.fa"
        File transcriptome_ref_fasta_index = "gs://broad-dsde-methods-long-reads-public/resources/gencode_v37/gencode.v37.pc_transcripts.fa.fai"
        File transcriptome_ref_fasta_dict = "gs://broad-dsde-methods-long-reads-public/resources/gencode_v37/gencode.v37.pc_transcripts.dict"

        # Default here is 0 because ccs uncorrected reads all seem to have RQ = -1.
        # All pathologically long reads also have RQ = -1.
        # This way we preserve the vast majority of the data, even if it has low quality.
        # We can filter it out at later steps.
        Float min_read_quality = 0.0
        Int max_reclamation_length = 60000

        String mas_seq_model = "mas15"


        # Set up some meta parameters here so we can adjust for when we want things to go VERY fast:
        Int primary_scatter_width = 50
        Int secondary_scatter_width = 10
    }

    parameter_meta {
        bam_file : "Input BAM file to process.  MUST already have been annotated by Longbow."
        sample_name : "The name of the sample to associate with the data in this workflow."
        gcs_out_root_dir : "Root output GCS folder in which to place results of this workflow."

        ref_fasta : "FASTA file containing the reference sequence to which the input data should be aligned before splitting into array elements."
        ref_fasta_index : "FASTA index file for the given ref_fasta file."
        ref_fasta_dict : "Sequence dictionary file for the given ref_fasta file."

        transcriptome_ref_fasta : "FASTA file containing the reference sequence to which the array elements should be aligned."
        transcriptome_ref_fasta_index : "FASTA index file for the given transcriptome_ref_fasta file."
        transcriptome_ref_fasta_dict : "Sequence dictionary file for the given transcriptome_ref_fasta file."

        min_read_quality : "[optional] Minimum read quality for reads to have to be included in our data (Default: 0.0)."
        max_reclamation_length : "[optional] Maximum length (in bases) that a read can be to attempt to reclaim from CCS rejection (Default: 60000)."

        primary_scatter_width : "[optional] Width to use for the primary scatter operation on this dataset (default: 50)."
        secondary_scatter_width : "[optional] Width to use for the secondary (nested) scatter operation on this dataset (default: 10)."
    }

    # Version of this workflow.
    String VERSION = "0.2"

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_01_WdlExecutionStartTimestamp { input: }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    # Define some attributes for later:
    RuntimeAttr disable_preemption_runtime_attrs = object {
        preemptible_tries: 0
    }
    RuntimeAttr filterReadsAttrs = object {
        cpu_cores: 4,
        preemptible_tries: 0
    }

    call PB.PBIndex as t_02_PbIndexInputBam {
        input:
            bam = bam_file
    }

    call PB.ShardLongReads as t_03_ShardLongReads {
        input:
            unaligned_bam = bam_file,
            unaligned_pbi = t_02_PbIndexInputBam.pbindex,
            prefix = sample_name + "_shard",
            num_shards = primary_scatter_width,
    }

    scatter (sharded_reads in t_03_ShardLongReads.unmapped_shards) {

        #######################################################################
        ## Setup
        #######################################################################

        String fbmrq_prefix = basename(sharded_reads, ".bam")

        #######################################################################
        ## Handle CCS Corrected Reads
        #######################################################################

        # 1 - filter the reads by the minimum read quality:
        call Utils.Bamtools as t_04_FilterByMinQual {
            input:
                bamfile = sharded_reads,
                prefix = fbmrq_prefix + "_good_reads",
                cmd = "filter",
                args = '-tag "rq":">=' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # 6: Longbow filter ccs reclaimable reads
        call LONGBOW.Filter as t_05_FilterCCSReads {
            input:
                bam = t_04_FilterByMinQual.bam_out,
                prefix = sample_name + "_ccs_corrected_subshard",
                model = mas_seq_model
        }

        call PB.PBIndex as t_06_PbIndexLongbowAnnotatedCCSPassedReads {
            input:
                bam = t_05_FilterCCSReads.passed_reads
        }

        # Shard these reads even wider so we can make sure we don't run out of memory:
        call PB.ShardLongReads as t_07_ShardCorrectedReads {
            input:
                unaligned_bam = t_05_FilterCCSReads.passed_reads,
                unaligned_pbi = t_06_PbIndexLongbowAnnotatedCCSPassedReads.pbindex,
                prefix = sample_name + "_ccs_corrected_longbow_annotated_subshard",
                num_shards = secondary_scatter_width,
        }

        # Segment our arrays into individual array elements:
        scatter (corrected_shard in t_07_ShardCorrectedReads.unmapped_shards) {
            call LONGBOW.Segment as t_08_SegmentCCSAnnotatedReads {
                input:
                    annotated_reads = corrected_shard,
                    model = mas_seq_model
            }
        }

        # Merge all outputs of Longbow Annotate / Segment:
        call Utils.MergeBams as t_09_MergeCCSArrayElements_1 {
            input:
                bams = t_08_SegmentCCSAnnotatedReads.segmented_bam,
                prefix = sample_name + "_ccs_corrected_ArrayElements_intermediate_1"
        }

        # Align our array elements to the transcriptome:
        call AR.Minimap2 as t_10_AlignCCSArrayElements {
            input:
                reads      = [ t_09_MergeCCSArrayElements_1.merged_bam ],
                ref_fasta  = transcriptome_ref_fasta,
                map_preset = "asm20"
        }

        # Align our array elements to the genome in splice-aware mode:
        call AR.Minimap2 as t_11_AlignCCSArrayElementsToGenome {
            input:
                reads      = [ t_09_MergeCCSArrayElements_1.merged_bam ],
                ref_fasta  = ref_fasta,
                map_preset = "splice:hq"
        }

        # We need to restore the annotations we created with the 10x tool to the transcriptome aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_12_RestoreAnnotationsToTranscriptomeAlignedCCSBam {
            input:
                annotated_bam_file = t_09_MergeCCSArrayElements_1.merged_bam,
                aligned_bam_file = t_10_AlignCCSArrayElements.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # We need to restore the annotations we created with the 10x tool to the genome aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_13_RestoreAnnotationsToGenomeAlignedCCSBam {
            input:
                annotated_bam_file = t_09_MergeCCSArrayElements_1.merged_bam,
                aligned_bam_file = t_11_AlignCCSArrayElementsToGenome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # To properly count our transcripts we must throw away the non-primary and unaligned reads:
        # Filter out all non-primary transcriptome alignments:
        call Utils.FilterReadsBySamFlags as t_14_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads {
            input:
                bam = t_12_RestoreAnnotationsToTranscriptomeAlignedCCSBam.output_bam,
                sam_flags = "2308",
                prefix = sample_name + "_ccs_corrected_ArrayElements_Annotated_Transcriptome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }
        # Filter out all non-primary genome alignments:
        call Utils.FilterReadsBySamFlags as t_15_RemoveUnmappedAndNonPrimaryCCSGenomeReads {
            input:
                bam = t_13_RestoreAnnotationsToGenomeAlignedCCSBam.output_bam,
                sam_flags = "2308",
                prefix = sample_name + "_ccs_corrected_ArrayElements_Annotated_Genome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }

        #######################################################################
        ## Handle CCS Uncorrected Reads
        #######################################################################

        # 1.5 - Get the "rejected" reads:
        call Utils.Bamtools as t_16_GetCcsRejectedReads {
            input:
                bamfile = sharded_reads,
                prefix = fbmrq_prefix + "_rejected_reads",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # 2 - Get reads we can reclaim:
        call Utils.Bamtools as t_17_ExtractCcsReclaimableReads {
            input:
                bamfile = sharded_reads,
                prefix = fbmrq_prefix + "_reads_for_ccs_reclamation",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '" -length "<=' + max_reclamation_length + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        call PB.PBIndex as t_18_PbIndexLongbowAnnotatedReclaimedReads {
            input:
                bam = t_17_ExtractCcsReclaimableReads.bam_out
        }

        # 6: Longbow filter ccs reclaimable reads
        call LONGBOW.Filter as t_19_FilterReclaimableReads {
            input:
                bam = t_17_ExtractCcsReclaimableReads.bam_out,
                bam_pbi = t_18_PbIndexLongbowAnnotatedReclaimedReads.pbindex,
                prefix = sample_name + "_subshard",
                model = mas_seq_model
        }

        call PB.PBIndex as t_20_PbIndexLongbowAnnotatedReclaimedPassedReads {
            input:
                bam = t_19_FilterReclaimableReads.passed_reads
        }

        # Shard these reads even wider so we can make sure we don't run out of memory:
        call PB.ShardLongReads as t_21_ShardReclaimedReads {
            input:
                unaligned_bam = t_19_FilterReclaimableReads.passed_reads,
                unaligned_pbi = t_20_PbIndexLongbowAnnotatedReclaimedPassedReads.pbindex,
                prefix = sample_name + "_ccs_reclaimed_longbow_annotated_subshard",
                num_shards = secondary_scatter_width,
        }

        # Segment our arrays into individual array elements:
        scatter (corrected_shard in t_21_ShardReclaimedReads.unmapped_shards) {
            call LONGBOW.Segment as t_22_SegmentReclaimedAnnotatedReads {
                input:
                    annotated_reads = corrected_shard,
                    model = mas_seq_model
            }
        }

        # Merge all outputs of Longbow Annotate / Segment:
        call Utils.MergeBams as t_23_MergeReclaimedArrayElements_1 {
            input:
                bams = t_22_SegmentReclaimedAnnotatedReads.segmented_bam,
                prefix = sample_name + "_ccs_reclaimed_ArrayElements_intermediate_1"
        }

        # Align our array elements to the transcriptome:
        call AR.Minimap2 as t_24_AlignReclaimedArrayElements {
            input:
                reads      = [ t_23_MergeReclaimedArrayElements_1.merged_bam ],
                ref_fasta  = transcriptome_ref_fasta,
                map_preset = "asm20"
        }

        # Align our array elements to the genome in splice-aware mode:
        call AR.Minimap2 as t_25_AlignReclaimedArrayElementsToGenome {
            input:
                reads      = [ t_23_MergeReclaimedArrayElements_1.merged_bam ],
                ref_fasta  = ref_fasta,
                map_preset = "splice:hq"
        }

        # We need to restore the annotations we created with the 10x tool to the transcriptome aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_26_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam {
            input:
                annotated_bam_file = t_23_MergeReclaimedArrayElements_1.merged_bam,
                aligned_bam_file = t_24_AlignReclaimedArrayElements.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # We need to restore the annotations we created with the 10x tool to the genome aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_27_RestoreAnnotationsToGenomeAlignedReclaimedBam {
            input:
                annotated_bam_file = t_23_MergeReclaimedArrayElements_1.merged_bam,
                aligned_bam_file = t_25_AlignReclaimedArrayElementsToGenome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # To properly count our transcripts we must throw away the non-primary and unaligned reads:
        # Filter out all non-primary transcriptome alignments:
        call Utils.FilterReadsBySamFlags as t_28_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads {
            input:
                bam = t_26_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam.output_bam,
                sam_flags = "2308",
                prefix = sample_name + "_ccs_reclaimed_ArrayElements_Annotated_Transcriptome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }
        # Filter out all non-primary genome alignments:
        call Utils.FilterReadsBySamFlags as t_29_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads {
            input:
                bam = t_27_RestoreAnnotationsToGenomeAlignedReclaimedBam.output_bam,
                sam_flags = "2308",
                prefix = sample_name + "_ccs_reclaimed_ArrayElements_Annotated_Genome_Aligned_PrimaryOnly",
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
    call Utils.MergeBams as t_30_MergeCCSReads { input: bams = t_04_FilterByMinQual.bam_out, prefix = sample_name + "_ccs_reads" }

    # Merge all CCS bams together for this Subread BAM:
    call Utils.MergeBams as t_31_MergePassedCCSReads { input: bams = t_05_FilterCCSReads.passed_reads, prefix = sample_name + "_ccs_longbow_passed_reads" }
    call Utils.MergeBams as t_32_MergeFailedCCSReads { input: bams = t_05_FilterCCSReads.failed_reads, prefix = sample_name + "_ccs_longbow_failed_reads" }
    call Utils.MergeBams as t_33_MergeCCSArrayElements { input: bams = t_09_MergeCCSArrayElements_1.merged_bam, prefix = sample_name + "_ccs_array_elements" }
    call Utils.MergeBams as t_34_MergeTranscriptomeAlignedCCSArrayElements { input: bams = t_12_RestoreAnnotationsToTranscriptomeAlignedCCSBam.output_bam, prefix = sample_name + "_ccs_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_35_MergeGenomeAlignedCCSArrayElements { input: bams = t_13_RestoreAnnotationsToGenomeAlignedCCSBam.output_bam, prefix = sample_name + "_ccs_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_36_MergePrimaryTranscriptomeAlignedCCSArrayElements { input: bams = t_14_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads.output_bam, prefix = sample_name + "_ccs_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_37_MergePrimaryGenomeAlignedCCSArrayElements { input: bams = t_15_RemoveUnmappedAndNonPrimaryCCSGenomeReads.output_bam, prefix = sample_name + "_ccs_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    #################################################
    # CCS Reclaimed:
    call Utils.MergeBams as t_38_MergeCCSReclaimableReads { input: bams = t_17_ExtractCcsReclaimableReads.bam_out, prefix = sample_name + "_reclaimable_reads" }

    # Merge all CCS bams together for this Subread BAM:
    call Utils.MergeBams as t_39_MergePassedCCSReclaimedReads { input: bams = t_19_FilterReclaimableReads.passed_reads, prefix = sample_name + "_reclaimable_longbow_passed_reads" }
    call Utils.MergeBams as t_40_MergeFailedCCSReclaimableReads { input: bams = t_19_FilterReclaimableReads.failed_reads, prefix = sample_name + "_reclaimed_longbow_failed_reads" }
    call Utils.MergeBams as t_41_MergeCCSReclaimedArrayElements { input: bams = t_23_MergeReclaimedArrayElements_1.merged_bam, prefix = sample_name + "_reclaimed_array_elements" }
    call Utils.MergeBams as t_42_MergeTranscriptomeAlignedCCSReclaimedArrayElements { input: bams = t_26_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam.output_bam, prefix = sample_name + "_reclaimed_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_43_MergeGenomeAlignedCCSReclaimedArrayElements { input: bams = t_27_RestoreAnnotationsToGenomeAlignedReclaimedBam.output_bam, prefix = sample_name + "_reclaimed_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_44_MergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements { input: bams = t_28_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads.output_bam, prefix = sample_name + "_reclaimed_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_45_MergePrimaryGenomeAlignedCCSReclaimedArrayElements { input: bams = t_29_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads.output_bam, prefix = sample_name + "_reclaimed_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    #################################################
    # Overall:
    call Utils.MergeBams as t_46_MergeAllLongbowPassedReads { input: bams = flatten([t_05_FilterCCSReads.passed_reads, t_19_FilterReclaimableReads.passed_reads]), prefix = sample_name + "_all_longbow_passed_reads" }
    call Utils.MergeBams as t_47_MergeAllLongbowFailedReads { input: bams = flatten([t_05_FilterCCSReads.failed_reads, t_19_FilterReclaimableReads.failed_reads]), prefix = sample_name + "_all_longbow_failed_reads" }

    call Utils.MergeBams as t_48_MergeAllArrayElements { input: bams = flatten([t_09_MergeCCSArrayElements_1.merged_bam, t_23_MergeReclaimedArrayElements_1.merged_bam]), prefix = sample_name + "_all_array_elements", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_49_MergeAllTranscriptomeAlignedArrayElements { input: bams = flatten([t_12_RestoreAnnotationsToTranscriptomeAlignedCCSBam.output_bam, t_26_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam.output_bam]), prefix = sample_name + "_all_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_50_MergeAllGenomeAlignedArrayElements { input: bams = flatten([t_13_RestoreAnnotationsToGenomeAlignedCCSBam.output_bam, t_27_RestoreAnnotationsToGenomeAlignedReclaimedBam.output_bam]), prefix = sample_name + "_all_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_51_MergeAllPrimaryTranscriptomeAlignedArrayElements { input: bams = flatten([t_14_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads.output_bam, t_28_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads.output_bam]), prefix = sample_name + "_all_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_52_MergeAllPrimaryGenomeAlignedArrayElements { input: bams = flatten([t_15_RemoveUnmappedAndNonPrimaryCCSGenomeReads.output_bam, t_29_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads.output_bam]), prefix = sample_name + "_all_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    #################################################
    #   ___      ____
    #  / _ \    / ___|
    # | | | |  | |
    # | |_| |  | |___
    #  \__\_\   \____|
    #
    #################################################

    # Get stats on CCS reads:
    call LONGBOW.Stats as t_53_CCS_longbow_stats {
        input:
            reads = t_30_MergeCCSReads.merged_bam,
            model = mas_seq_model,
            prefix = sample_name + "_CCS_Corrected"
    }

    # Get stats on Reclaimable reads:
    call LONGBOW.Stats as t_54_Reclaimable_longbow_stats {
        input:
            reads = t_38_MergeCCSReclaimableReads.merged_bam,
            model = mas_seq_model,
            prefix = sample_name + "_CCS_Reclaimable"
    }

    # Get stats on Reclaimed reads:
    call LONGBOW.Stats as t_55_Reclaimed_longbow_stats {
        input:
            reads = t_39_MergePassedCCSReclaimedReads.merged_bam,
            model = mas_seq_model,
            prefix = sample_name + "_CCS_Reclaimed"
    }

    # Get stats on All reads (overall stats):
    call LONGBOW.Stats as t_56_Passed_longbow_stats {
        input:
            reads = t_46_MergeAllLongbowPassedReads.merged_bam,
            model = mas_seq_model,
            prefix = sample_name + "_All_Longbow_Passed"
    }

    # Get stats on All reads (overall stats):
    call LONGBOW.Stats as t_57_Overall_longbow_stats {
        input:
            reads = bam_file,
            model = mas_seq_model,
            prefix = sample_name + "_Overall"
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
    String base_out_dir = outdir + "/" + sample_name+ "/" + t_01_WdlExecutionStartTimestamp.timestamp_string

    String array_element_dir = base_out_dir + "/array_elements"
    String arrays_dir = base_out_dir + "/array_reads"

    ##############################################################################################################
    # Finalize the CCS Array Reads
    call FF.FinalizeToDir as t_58_FinalizeCCSReads {
        input:
            files = [
                t_30_MergeCCSReads.merged_bam,
                t_30_MergeCCSReads.merged_bai,
            ],
            outdir = arrays_dir,
            keyfile = t_57_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Reclaimed Array Reads
    call FF.FinalizeToDir as t_59_FinalizeCCSReclaimedReads {
        input:
            files = [
                t_38_MergeCCSReclaimableReads.merged_bam,
                t_38_MergeCCSReclaimableReads.merged_bai,
                t_39_MergePassedCCSReclaimedReads.merged_bam,
                t_39_MergePassedCCSReclaimedReads.merged_bai,
                t_40_MergeFailedCCSReclaimableReads.merged_bam,
                t_40_MergeFailedCCSReclaimableReads.merged_bai,
            ],
            outdir = arrays_dir,
            keyfile = t_57_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Overall Array Reads
    call FF.FinalizeToDir as t_60_FinalizeOverallCombinedReads {
        input:
            files = [
                t_46_MergeAllLongbowPassedReads.merged_bam,
                t_46_MergeAllLongbowPassedReads.merged_bai,
                t_47_MergeAllLongbowFailedReads.merged_bam,
                t_47_MergeAllLongbowFailedReads.merged_bai,
            ],
            outdir = arrays_dir,
            keyfile = t_57_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the CCS Array Element files
    call FF.FinalizeToDir as t_61_FinalizeCCSArrayElements {
        input:
            files = [
                t_33_MergeCCSArrayElements.merged_bam,
                t_33_MergeCCSArrayElements.merged_bai,
                t_34_MergeTranscriptomeAlignedCCSArrayElements.merged_bam,
                t_34_MergeTranscriptomeAlignedCCSArrayElements.merged_bai,
                t_35_MergeGenomeAlignedCCSArrayElements.merged_bam,
                t_35_MergeGenomeAlignedCCSArrayElements.merged_bai,
                t_36_MergePrimaryTranscriptomeAlignedCCSArrayElements.merged_bam,
                t_36_MergePrimaryTranscriptomeAlignedCCSArrayElements.merged_bai,
                t_37_MergePrimaryGenomeAlignedCCSArrayElements.merged_bam,
                t_37_MergePrimaryGenomeAlignedCCSArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = t_57_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Reclaimed Array Element files
    call FF.FinalizeToDir as t_62_FinalizeCCSRecaimedArrayElements {
        input:
            files = [
                t_41_MergeCCSReclaimedArrayElements.merged_bam,
                t_41_MergeCCSReclaimedArrayElements.merged_bai,
                t_42_MergeTranscriptomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_42_MergeTranscriptomeAlignedCCSReclaimedArrayElements.merged_bai,
                t_43_MergeGenomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_43_MergeGenomeAlignedCCSReclaimedArrayElements.merged_bai,
                t_44_MergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_44_MergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements.merged_bai,
                t_45_MergePrimaryGenomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_45_MergePrimaryGenomeAlignedCCSReclaimedArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = t_57_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Overall Array Element files
    call FF.FinalizeToDir as t_63_FinalizeOverallArrayElements {
        input:
            files = [
                t_48_MergeAllArrayElements.merged_bam,
                t_48_MergeAllArrayElements.merged_bai,
                t_49_MergeAllTranscriptomeAlignedArrayElements.merged_bam,
                t_49_MergeAllTranscriptomeAlignedArrayElements.merged_bai,
                t_50_MergeAllGenomeAlignedArrayElements.merged_bam,
                t_50_MergeAllGenomeAlignedArrayElements.merged_bai,
                t_51_MergeAllPrimaryTranscriptomeAlignedArrayElements.merged_bam,
                t_51_MergeAllPrimaryTranscriptomeAlignedArrayElements.merged_bai,
                t_52_MergeAllPrimaryGenomeAlignedArrayElements.merged_bam,
                t_52_MergeAllPrimaryGenomeAlignedArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = t_57_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the different stats files to separate directories
    call FF.FinalizeToDir as t_64_FinalizeCCSLongbowStats {
        input:
            files = [
                t_53_CCS_longbow_stats.summary_stats,
                t_53_CCS_longbow_stats.array_length_counts_plot_png,
                t_53_CCS_longbow_stats.array_length_counts_plot_svg,
                t_53_CCS_longbow_stats.ligation_heatmap_nn_png,
                t_53_CCS_longbow_stats.ligation_heatmap_nn_svg,
                t_53_CCS_longbow_stats.ligation_heatmap_png,
                t_53_CCS_longbow_stats.ligation_heatmap_svg,
                t_53_CCS_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_53_CCS_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_53_CCS_longbow_stats.ligation_heatmap_reduced_png,
                t_53_CCS_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/CCS_Corrected/",
            keyfile = t_57_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_65_FinalizeReclaimableLongbowStats {
        input:
            files = [
                t_54_Reclaimable_longbow_stats.summary_stats,
                t_54_Reclaimable_longbow_stats.array_length_counts_plot_png,
                t_54_Reclaimable_longbow_stats.array_length_counts_plot_svg,
                t_54_Reclaimable_longbow_stats.ligation_heatmap_nn_png,
                t_54_Reclaimable_longbow_stats.ligation_heatmap_nn_svg,
                t_54_Reclaimable_longbow_stats.ligation_heatmap_png,
                t_54_Reclaimable_longbow_stats.ligation_heatmap_svg,
                t_54_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_54_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_54_Reclaimable_longbow_stats.ligation_heatmap_reduced_png,
                t_54_Reclaimable_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/CCS_Reclaimable/",
            keyfile = t_57_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_66_FinalizeReclaimedLongbowStats {
        input:
            files = [
                t_55_Reclaimed_longbow_stats.summary_stats,
                t_55_Reclaimed_longbow_stats.array_length_counts_plot_png,
                t_55_Reclaimed_longbow_stats.array_length_counts_plot_svg,
                t_55_Reclaimed_longbow_stats.ligation_heatmap_nn_png,
                t_55_Reclaimed_longbow_stats.ligation_heatmap_nn_svg,
                t_55_Reclaimed_longbow_stats.ligation_heatmap_png,
                t_55_Reclaimed_longbow_stats.ligation_heatmap_svg,
                t_55_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_55_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_55_Reclaimed_longbow_stats.ligation_heatmap_reduced_png,
                t_55_Reclaimed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/CCS_Reclaimed/",
            keyfile = t_57_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_67_FinalizeOverallLongbowStats {
        input:
            files = [
                t_57_Overall_longbow_stats.summary_stats,
                t_57_Overall_longbow_stats.array_length_counts_plot_png,
                t_57_Overall_longbow_stats.array_length_counts_plot_svg,
                t_57_Overall_longbow_stats.ligation_heatmap_nn_png,
                t_57_Overall_longbow_stats.ligation_heatmap_nn_svg,
                t_57_Overall_longbow_stats.ligation_heatmap_png,
                t_57_Overall_longbow_stats.ligation_heatmap_svg,
                t_57_Overall_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_57_Overall_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_57_Overall_longbow_stats.ligation_heatmap_reduced_png,
                t_57_Overall_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/Overall/",
            keyfile = t_57_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_68_FinalizeAllPassedLongbowStats {
        input:
            files = [
                t_56_Passed_longbow_stats.summary_stats,
                t_56_Passed_longbow_stats.array_length_counts_plot_png,
                t_56_Passed_longbow_stats.array_length_counts_plot_svg,
                t_56_Passed_longbow_stats.ligation_heatmap_nn_png,
                t_56_Passed_longbow_stats.ligation_heatmap_nn_svg,
                t_56_Passed_longbow_stats.ligation_heatmap_png,
                t_56_Passed_longbow_stats.ligation_heatmap_svg,
                t_56_Passed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_56_Passed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_56_Passed_longbow_stats.ligation_heatmap_reduced_png,
                t_56_Passed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/All_Longbow_Passed/",
            keyfile = t_57_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good
    call FF.WriteCompletionFile as t_69_WriteCompletionFile {
        input:
            outdir = base_out_dir + "/",
            keyfile =  t_57_Overall_longbow_stats.summary_stats
    }
}
