version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF
import "tasks/AlignReads.wdl" as AR
import "tasks/Ten_X_Tool.wdl" as TENX
import "tasks/Longbow.wdl" as LONGBOW
import "tasks/Structs.wdl"

workflow MasSeqAlignmentOnlyNoAnnotate {

    meta {
        description : "This workflow will process and align MAS-seq data.  No UMI / CBC calculations are performed."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File reads_bam
        String sample_name
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/MasSeqAlignmentOnly"

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
        Int primary_scatter_width = 200
        Int secondary_scatter_width = 10

        String longbow_docker_version = "us.gcr.io/broad-dsp-lrma/lr-longbow:0.4.3"
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

    RuntimeAttr new_longbow_attrs = object {
        docker: longbow_docker_version
    }

    File read_pbi = sub(reads_bam, ".bam$", ".bam.pbi")
    call PB.ShardLongReads as t_02_ShardLongReads {
        input:
            unaligned_bam = reads_bam,
            unaligned_pbi = read_pbi,
            prefix = sample_name + "_shard",
            num_shards = primary_scatter_width,
    }
    
    scatter (sharded_reads in t_02_ShardLongReads.unmapped_shards) {

        #######################################################################
        ## Setup
        #######################################################################

        String fbmrq_prefix = basename(sharded_reads, ".bam")

        #######################################################################
        ## Handle CCS Corrected Reads
        #######################################################################

        # 1 - filter the reads by the minimum read quality:
        call Utils.Bamtools as t_03_FilterByMinQual {
            input:
                bamfile = sharded_reads,
                prefix = fbmrq_prefix + "_good_reads",
                cmd = "filter",
                args = '-tag "rq":">=' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # 6: Longbow filter ccs reclaimable reads
        call LONGBOW.Filter as t_04_FilterCCSReads {
            input:
                bam = t_03_FilterByMinQual.bam_out,
                prefix = sample_name + "_ccs_corrected_subshard",
                model = mas_seq_model,
                runtime_attr_override = new_longbow_attrs,
        }

        call PB.PBIndex as t_05_PbIndexLongbowAnnotatedCCSPassedReads {
            input:
                bam = t_04_FilterCCSReads.passed_reads
        }

        # Shard these reads even wider so we can make sure we don't run out of memory:
        call PB.ShardLongReads as t_06_ShardCorrectedReads {
            input:
                unaligned_bam = t_04_FilterCCSReads.passed_reads,
                unaligned_pbi = t_05_PbIndexLongbowAnnotatedCCSPassedReads.pbindex,
                prefix = sample_name + "_ccs_corrected_longbow_annotated_subshard",
                num_shards = secondary_scatter_width,
        }

        # Segment our arrays into individual array elements:
        scatter (corrected_shard in t_06_ShardCorrectedReads.unmapped_shards) {
            call LONGBOW.Segment as t_07_SegmentCCSAnnotatedReads {
                input:
                    annotated_reads = corrected_shard,
                    model = mas_seq_model,
                    runtime_attr_override = new_longbow_attrs,
            }
        }

        # Merge all outputs of Longbow Annotate / Segment:
        call Utils.MergeBams as t_08_MergeCCSArrayElements_1 {
            input:
                bams = t_07_SegmentCCSAnnotatedReads.segmented_bam,
                prefix = sample_name + "_ccs_corrected_ArrayElements_intermediate_1"
        }

        # Align our array elements to the transcriptome:
        call AR.Minimap2 as t_09_AlignCCSArrayElements {
            input:
                reads      = [ t_08_MergeCCSArrayElements_1.merged_bam ],
                ref_fasta  = transcriptome_ref_fasta,
                map_preset = "asm20"
        }

        # Align our array elements to the genome in splice-aware mode:
        call AR.Minimap2 as t_10_AlignCCSArrayElementsToGenome {
            input:
                reads      = [ t_08_MergeCCSArrayElements_1.merged_bam ],
                ref_fasta  = ref_fasta,
                map_preset = "splice:hq"
        }

        # We need to restore the annotations we created with the 10x tool to the transcriptome aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_11_RestoreAnnotationsToTranscriptomeAlignedCCSBam {
            input:
                annotated_bam_file = t_08_MergeCCSArrayElements_1.merged_bam,
                aligned_bam_file = t_09_AlignCCSArrayElements.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # We need to restore the annotations we created with the 10x tool to the genome aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_12_RestoreAnnotationsToGenomeAlignedCCSBam {
            input:
                annotated_bam_file = t_08_MergeCCSArrayElements_1.merged_bam,
                aligned_bam_file = t_10_AlignCCSArrayElementsToGenome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # To properly count our transcripts we must throw away the non-primary and unaligned reads:
        # Filter out all non-primary transcriptome alignments:
        call Utils.FilterReadsBySamFlags as t_13_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads {
            input:
                bam = t_11_RestoreAnnotationsToTranscriptomeAlignedCCSBam.output_bam,
                sam_flags = "2308",
                prefix = sample_name + "_ccs_corrected_ArrayElements_Annotated_Transcriptome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }
        # Filter out all non-primary genome alignments:
        call Utils.FilterReadsBySamFlags as t_14_RemoveUnmappedAndNonPrimaryCCSGenomeReads {
            input:
                bam = t_12_RestoreAnnotationsToGenomeAlignedCCSBam.output_bam,
                sam_flags = "2308",
                prefix = sample_name + "_ccs_corrected_ArrayElements_Annotated_Genome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }

        #######################################################################
        ## Handle CCS Uncorrected Reads
        #######################################################################

        # 2 - Get reads we can reclaim:
        call Utils.Bamtools as t_15_ExtractCcsReclaimableReads {
            input:
                bamfile = sharded_reads,
                prefix = fbmrq_prefix + "_reads_for_ccs_reclamation",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '" -length "<=' + max_reclamation_length + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        call PB.PBIndex as t_16_PbIndexLongbowAnnotatedReclaimedReads {
            input:
                bam = t_15_ExtractCcsReclaimableReads.bam_out,
        }

        # 6: Longbow filter ccs reclaimable reads
        call LONGBOW.Filter as t_17_FilterReclaimableReads {
            input:
                bam = t_15_ExtractCcsReclaimableReads.bam_out,
                bam_pbi = t_16_PbIndexLongbowAnnotatedReclaimedReads.pbindex,
                prefix = sample_name + "_subshard",
                model = mas_seq_model,
                runtime_attr_override = new_longbow_attrs,
        }

        call PB.PBIndex as t_18_PbIndexLongbowAnnotatedReclaimedPassedReads {
            input:
                bam = t_17_FilterReclaimableReads.passed_reads
        }

        # Shard these reads even wider so we can make sure we don't run out of memory:
        call PB.ShardLongReads as t_19_ShardReclaimedReads {
            input:
                unaligned_bam = t_17_FilterReclaimableReads.passed_reads,
                unaligned_pbi = t_18_PbIndexLongbowAnnotatedReclaimedPassedReads.pbindex,
                prefix = sample_name + "_ccs_reclaimed_longbow_annotated_subshard",
                num_shards = secondary_scatter_width,
        }

        # Segment our arrays into individual array elements:
        scatter (corrected_shard in t_19_ShardReclaimedReads.unmapped_shards) {
            call LONGBOW.Segment as t_20_SegmentReclaimedAnnotatedReads {
                input:
                    annotated_reads = corrected_shard,
                    model = mas_seq_model,
                    runtime_attr_override = new_longbow_attrs,
            }
        }

        # Merge all outputs of Longbow Annotate / Segment:
        call Utils.MergeBams as t_21_MergeReclaimedArrayElements_1 {
            input:
                bams = t_20_SegmentReclaimedAnnotatedReads.segmented_bam,
                prefix = sample_name + "_ccs_reclaimed_ArrayElements_intermediate_1"
        }

        # Align our array elements to the transcriptome:
        call AR.Minimap2 as t_22_AlignReclaimedArrayElements {
            input:
                reads      = [ t_21_MergeReclaimedArrayElements_1.merged_bam ],
                ref_fasta  = transcriptome_ref_fasta,
                map_preset = "asm20"
        }

        # Align our array elements to the genome in splice-aware mode:
        call AR.Minimap2 as t_23_AlignReclaimedArrayElementsToGenome {
            input:
                reads      = [ t_21_MergeReclaimedArrayElements_1.merged_bam ],
                ref_fasta  = ref_fasta,
                map_preset = "splice:hq"
        }

        # We need to restore the annotations we created with the 10x tool to the transcriptome aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_24_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam {
            input:
                annotated_bam_file = t_21_MergeReclaimedArrayElements_1.merged_bam,
                aligned_bam_file = t_22_AlignReclaimedArrayElements.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # We need to restore the annotations we created with the 10x tool to the genome aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_25_RestoreAnnotationsToGenomeAlignedReclaimedBam {
            input:
                annotated_bam_file = t_21_MergeReclaimedArrayElements_1.merged_bam,
                aligned_bam_file = t_23_AlignReclaimedArrayElementsToGenome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # To properly count our transcripts we must throw away the non-primary and unaligned reads:
        # Filter out all non-primary transcriptome alignments:
        call Utils.FilterReadsBySamFlags as t_26_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads {
            input:
                bam = t_24_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam.output_bam,
                sam_flags = "2308",
                prefix = sample_name + "_ccs_reclaimed_ArrayElements_Annotated_Transcriptome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }
        # Filter out all non-primary genome alignments:
        call Utils.FilterReadsBySamFlags as t_27_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads {
            input:
                bam = t_25_RestoreAnnotationsToGenomeAlignedReclaimedBam.output_bam,
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
    call Utils.MergeBams as t_28_MergeCCSReads { input: bams = t_03_FilterByMinQual.bam_out, prefix = sample_name + "_ccs_reads" }

    # Merge all CCS bams together for this Subread BAM:
    call Utils.MergeBams as t_29_MergePassedCCSReads { input: bams = t_04_FilterCCSReads.passed_reads, prefix = sample_name + "_ccs_longbow_passed_reads" }
    call Utils.MergeBams as t_30_MergeFailedCCSReads { input: bams = t_04_FilterCCSReads.failed_reads, prefix = sample_name + "_ccs_longbow_failed_reads" }
    call Utils.MergeBams as t_31_MergeCCSArrayElements { input: bams = t_08_MergeCCSArrayElements_1.merged_bam, prefix = sample_name + "_ccs_array_elements" }
    call Utils.MergeBams as t_32_MergeTranscriptomeAlignedCCSArrayElements { input: bams = t_11_RestoreAnnotationsToTranscriptomeAlignedCCSBam.output_bam, prefix = sample_name + "_ccs_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_33_MergeGenomeAlignedCCSArrayElements { input: bams = t_12_RestoreAnnotationsToGenomeAlignedCCSBam.output_bam, prefix = sample_name + "_ccs_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_34_MergePrimaryTranscriptomeAlignedCCSArrayElements { input: bams = t_13_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads.output_bam, prefix = sample_name + "_ccs_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_35_MergePrimaryGenomeAlignedCCSArrayElements { input: bams = t_14_RemoveUnmappedAndNonPrimaryCCSGenomeReads.output_bam, prefix = sample_name + "_ccs_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    #################################################
    # CCS Reclaimed:
    call Utils.MergeBams as t_36_MergeCCSReclaimableReads { input: bams = t_15_ExtractCcsReclaimableReads.bam_out, prefix = sample_name + "_reclaimable_reads" }

    # Merge all CCS bams together for this Subread BAM:
    call Utils.MergeBams as t_37_MergePassedCCSReclaimedReads { input: bams = t_17_FilterReclaimableReads.passed_reads, prefix = sample_name + "_reclaimable_longbow_passed_reads" }
    call Utils.MergeBams as t_38_MergeFailedCCSReclaimableReads { input: bams = t_17_FilterReclaimableReads.failed_reads, prefix = sample_name + "_reclaimed_longbow_failed_reads" }
    call Utils.MergeBams as t_39_MergeCCSReclaimedArrayElements { input: bams = t_21_MergeReclaimedArrayElements_1.merged_bam, prefix = sample_name + "_reclaimed_array_elements" }
    call Utils.MergeBams as t_40_MergeTranscriptomeAlignedCCSReclaimedArrayElements { input: bams = t_24_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam.output_bam, prefix = sample_name + "_reclaimed_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_41_MergeGenomeAlignedCCSReclaimedArrayElements { input: bams = t_25_RestoreAnnotationsToGenomeAlignedReclaimedBam.output_bam, prefix = sample_name + "_reclaimed_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_42_MergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements { input: bams = t_26_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads.output_bam, prefix = sample_name + "_reclaimed_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_43_MergePrimaryGenomeAlignedCCSReclaimedArrayElements { input: bams = t_27_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads.output_bam, prefix = sample_name + "_reclaimed_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    #################################################
    # Overall:
    call Utils.MergeBams as t_44_MergeAllLongbowPassedReads { input: bams = flatten([t_04_FilterCCSReads.passed_reads, t_17_FilterReclaimableReads.passed_reads]), prefix = sample_name + "_all_longbow_passed_reads" }
    call Utils.MergeBams as t_45_MergeAllLongbowFailedReads { input: bams = flatten([t_04_FilterCCSReads.failed_reads, t_17_FilterReclaimableReads.failed_reads]), prefix = sample_name + "_all_longbow_failed_reads" }

    call Utils.MergeBams as t_46_MergeAllArrayElements { input: bams = flatten([t_08_MergeCCSArrayElements_1.merged_bam, t_21_MergeReclaimedArrayElements_1.merged_bam]), prefix = sample_name + "_all_array_elements", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_47_MergeAllTranscriptomeAlignedArrayElements { input: bams = flatten([t_11_RestoreAnnotationsToTranscriptomeAlignedCCSBam.output_bam, t_24_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam.output_bam]), prefix = sample_name + "_all_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_48_MergeAllGenomeAlignedArrayElements { input: bams = flatten([t_12_RestoreAnnotationsToGenomeAlignedCCSBam.output_bam, t_25_RestoreAnnotationsToGenomeAlignedReclaimedBam.output_bam]), prefix = sample_name + "_all_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_49_MergeAllPrimaryTranscriptomeAlignedArrayElements { input: bams = flatten([t_13_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads.output_bam, t_26_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads.output_bam]), prefix = sample_name + "_all_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_50_MergeAllPrimaryGenomeAlignedArrayElements { input: bams = flatten([t_14_RemoveUnmappedAndNonPrimaryCCSGenomeReads.output_bam, t_27_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads.output_bam]), prefix = sample_name + "_all_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    #################################################
    #   ___      ____
    #  / _ \    / ___|
    # | | | |  | |
    # | |_| |  | |___
    #  \__\_\   \____|
    #
    #################################################

    # Get stats on CCS reads:
    call LONGBOW.Stats as t_51_CCS_longbow_stats {
        input:
            reads = t_28_MergeCCSReads.merged_bam,
            model = mas_seq_model,
            prefix = sample_name + "_CCS_Corrected",
            runtime_attr_override = new_longbow_attrs,
    }

    # Get stats on Reclaimable reads:
    call LONGBOW.Stats as t_52_Reclaimable_longbow_stats {
        input:
            reads = t_36_MergeCCSReclaimableReads.merged_bam,
            model = mas_seq_model,
            prefix = sample_name + "_CCS_Reclaimable",
            runtime_attr_override = new_longbow_attrs,
    }

    # Get stats on Reclaimed reads:
    call LONGBOW.Stats as t_53_Reclaimed_longbow_stats {
        input:
            reads = t_37_MergePassedCCSReclaimedReads.merged_bam,
            model = mas_seq_model,
            prefix = sample_name + "_CCS_Reclaimed",
            runtime_attr_override = new_longbow_attrs,
    }

    # Get stats on All reads (overall stats):
    call LONGBOW.Stats as t_54_Passed_longbow_stats {
        input:
            reads = t_44_MergeAllLongbowPassedReads.merged_bam,
            model = mas_seq_model,
            prefix = sample_name + "_All_Longbow_Passed",
            runtime_attr_override = new_longbow_attrs,
    }

    # Get stats on All reads (overall stats):
    call LONGBOW.Stats as t_55_Overall_longbow_stats {
        input:
            reads = reads_bam,
            model = mas_seq_model,
            prefix = sample_name + "_Overall",
            runtime_attr_override = new_longbow_attrs,
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
    String base_out_dir = outdir + "/" + sample_name + "/" + t_01_WdlExecutionStartTimestamp.timestamp_string

    String array_element_dir = base_out_dir + "/array_elements"
    String arrays_dir = base_out_dir + "/array_reads"

    ##############################################################################################################
    # Finalize the CCS Array Reads
    call FF.FinalizeToDir as t_56_FinalizeCCSReads {
        input:
            files = [
                t_28_MergeCCSReads.merged_bam,
                t_28_MergeCCSReads.merged_bai,
            ],
            outdir = arrays_dir,
            keyfile = t_55_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Reclaimed Array Reads
    call FF.FinalizeToDir as t_57_FinalizeCCSReclaimedReads {
        input:
            files = [
                t_36_MergeCCSReclaimableReads.merged_bam,
                t_36_MergeCCSReclaimableReads.merged_bai,
                t_37_MergePassedCCSReclaimedReads.merged_bam,
                t_37_MergePassedCCSReclaimedReads.merged_bai,
                t_38_MergeFailedCCSReclaimableReads.merged_bam,
                t_38_MergeFailedCCSReclaimableReads.merged_bai,
            ],
            outdir = arrays_dir,
            keyfile = t_55_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Overall Array Reads
    call FF.FinalizeToDir as t_58_FinalizeOverallCombinedReads {
        input:
            files = [
                t_44_MergeAllLongbowPassedReads.merged_bam,
                t_44_MergeAllLongbowPassedReads.merged_bai,
                t_45_MergeAllLongbowFailedReads.merged_bam,
                t_45_MergeAllLongbowFailedReads.merged_bai,
            ],
            outdir = arrays_dir,
            keyfile = t_55_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the CCS Array Element files
    call FF.FinalizeToDir as t_59_FinalizeCCSArrayElements {
        input:
            files = [
                t_31_MergeCCSArrayElements.merged_bam,
                t_31_MergeCCSArrayElements.merged_bai,
                t_32_MergeTranscriptomeAlignedCCSArrayElements.merged_bam,
                t_32_MergeTranscriptomeAlignedCCSArrayElements.merged_bai,
                t_33_MergeGenomeAlignedCCSArrayElements.merged_bam,
                t_33_MergeGenomeAlignedCCSArrayElements.merged_bai,
                t_34_MergePrimaryTranscriptomeAlignedCCSArrayElements.merged_bam,
                t_34_MergePrimaryTranscriptomeAlignedCCSArrayElements.merged_bai,
                t_35_MergePrimaryGenomeAlignedCCSArrayElements.merged_bam,
                t_35_MergePrimaryGenomeAlignedCCSArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = t_55_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Reclaimed Array Element files
    call FF.FinalizeToDir as t_60_FinalizeCCSRecaimedArrayElements {
        input:
            files = [
                t_39_MergeCCSReclaimedArrayElements.merged_bam,
                t_39_MergeCCSReclaimedArrayElements.merged_bai,
                t_40_MergeTranscriptomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_40_MergeTranscriptomeAlignedCCSReclaimedArrayElements.merged_bai,
                t_41_MergeGenomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_41_MergeGenomeAlignedCCSReclaimedArrayElements.merged_bai,
                t_42_MergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_42_MergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements.merged_bai,
                t_43_MergePrimaryGenomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_43_MergePrimaryGenomeAlignedCCSReclaimedArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = t_55_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Overall Array Element files
    call FF.FinalizeToDir as t_61_FinalizeOverallArrayElements {
        input:
            files = [
                t_46_MergeAllArrayElements.merged_bam,
                t_46_MergeAllArrayElements.merged_bai,
                t_47_MergeAllTranscriptomeAlignedArrayElements.merged_bam,
                t_47_MergeAllTranscriptomeAlignedArrayElements.merged_bai,
                t_48_MergeAllGenomeAlignedArrayElements.merged_bam,
                t_48_MergeAllGenomeAlignedArrayElements.merged_bai,
                t_49_MergeAllPrimaryTranscriptomeAlignedArrayElements.merged_bam,
                t_49_MergeAllPrimaryTranscriptomeAlignedArrayElements.merged_bai,
                t_50_MergeAllPrimaryGenomeAlignedArrayElements.merged_bam,
                t_50_MergeAllPrimaryGenomeAlignedArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = t_55_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the different stats files to separate directories
    call FF.FinalizeToDir as t_62_FinalizeCCSLongbowStats {
        input:
            files = [
                t_51_CCS_longbow_stats.summary_stats,
                t_51_CCS_longbow_stats.array_length_counts_plot_png,
                t_51_CCS_longbow_stats.array_length_counts_plot_svg,
                t_51_CCS_longbow_stats.ligation_heatmap_nn_png,
                t_51_CCS_longbow_stats.ligation_heatmap_nn_svg,
                t_51_CCS_longbow_stats.ligation_heatmap_png,
                t_51_CCS_longbow_stats.ligation_heatmap_svg,
                t_51_CCS_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_51_CCS_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_51_CCS_longbow_stats.ligation_heatmap_reduced_png,
                t_51_CCS_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/CCS_Corrected/",
            keyfile = t_55_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_63_FinalizeReclaimableLongbowStats {
        input:
            files = [
                t_52_Reclaimable_longbow_stats.summary_stats,
                t_52_Reclaimable_longbow_stats.array_length_counts_plot_png,
                t_52_Reclaimable_longbow_stats.array_length_counts_plot_svg,
                t_52_Reclaimable_longbow_stats.ligation_heatmap_nn_png,
                t_52_Reclaimable_longbow_stats.ligation_heatmap_nn_svg,
                t_52_Reclaimable_longbow_stats.ligation_heatmap_png,
                t_52_Reclaimable_longbow_stats.ligation_heatmap_svg,
                t_52_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_52_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_52_Reclaimable_longbow_stats.ligation_heatmap_reduced_png,
                t_52_Reclaimable_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/CCS_Reclaimable/",
            keyfile = t_55_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_64_FinalizeReclaimedLongbowStats {
        input:
            files = [
                t_53_Reclaimed_longbow_stats.summary_stats,
                t_53_Reclaimed_longbow_stats.array_length_counts_plot_png,
                t_53_Reclaimed_longbow_stats.array_length_counts_plot_svg,
                t_53_Reclaimed_longbow_stats.ligation_heatmap_nn_png,
                t_53_Reclaimed_longbow_stats.ligation_heatmap_nn_svg,
                t_53_Reclaimed_longbow_stats.ligation_heatmap_png,
                t_53_Reclaimed_longbow_stats.ligation_heatmap_svg,
                t_53_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_53_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_53_Reclaimed_longbow_stats.ligation_heatmap_reduced_png,
                t_53_Reclaimed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/CCS_Reclaimed/",
            keyfile = t_55_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_65_FinalizeOverallLongbowStats {
        input:
            files = [
                t_55_Overall_longbow_stats.summary_stats,
                t_55_Overall_longbow_stats.array_length_counts_plot_png,
                t_55_Overall_longbow_stats.array_length_counts_plot_svg,
                t_55_Overall_longbow_stats.ligation_heatmap_nn_png,
                t_55_Overall_longbow_stats.ligation_heatmap_nn_svg,
                t_55_Overall_longbow_stats.ligation_heatmap_png,
                t_55_Overall_longbow_stats.ligation_heatmap_svg,
                t_55_Overall_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_55_Overall_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_55_Overall_longbow_stats.ligation_heatmap_reduced_png,
                t_55_Overall_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/Overall/",
            keyfile = t_55_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_66_FinalizeAllPassedLongbowStats {
        input:
            files = [
                t_54_Passed_longbow_stats.summary_stats,
                t_54_Passed_longbow_stats.array_length_counts_plot_png,
                t_54_Passed_longbow_stats.array_length_counts_plot_svg,
                t_54_Passed_longbow_stats.ligation_heatmap_nn_png,
                t_54_Passed_longbow_stats.ligation_heatmap_nn_svg,
                t_54_Passed_longbow_stats.ligation_heatmap_png,
                t_54_Passed_longbow_stats.ligation_heatmap_svg,
                t_54_Passed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_54_Passed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_54_Passed_longbow_stats.ligation_heatmap_reduced_png,
                t_54_Passed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = base_out_dir + "/longbow_stats/All_Longbow_Passed/",
            keyfile = t_55_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good
    call FF.WriteCompletionFile as t_67_WriteCompletionFile {
        input:
            outdir = base_out_dir + "/",
            keyfile =  t_55_Overall_longbow_stats.summary_stats
    }
}
