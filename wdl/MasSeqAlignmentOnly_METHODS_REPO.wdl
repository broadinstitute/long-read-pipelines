version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    String? disk_type
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

workflow MasSeqAlignmentOnly {

    meta {
        description : "This workflow will process and align MAS-seq data.  No UMI / CBC calculations are performed."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String gcs_input_dir

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

        String? sample_name

        # Set up some meta parameters here so we can adjust for when we want things to go VERY fast:
        Int primary_scatter_width = 50
        Int secondary_scatter_width = 10
    }

    parameter_meta {
        gcs_input_dir : "Input folder on GCS in which to search for BAM files to process."

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

    # Version of this workflow.
    String VERSION = "0.2"

    call FindBams as t_01_FindBams { input: gcs_input_dir = gcs_input_dir }

    # Check here if we found ccs bams or subread bams:
    Boolean use_subreads = t_01_FindBams.has_subreads

    # Make sure we have **EXACTLY** one bam file to run on:
    if (length(t_01_FindBams.ccs_bams) != 1) {
        call FailWithWarning as t_02_FailOnMultiBamFiles { input: warning = "Error: Multiple BAM files found.  Cannot continue!" }
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
    File reads_bam = t_01_FindBams.ccs_bams[0]

    call GetRunInfo as t_03_GetRunInfo { input: subread_bam = reads_bam }

    String SM  = select_first([sample_name, t_03_GetRunInfo.run_info["SM"]])
    String PL  = "PACBIO"
    String PU  = t_03_GetRunInfo.run_info["PU"]
    String DT  = t_03_GetRunInfo.run_info["DT"]
    String ID  = PU
    String DS  = t_03_GetRunInfo.run_info["DS"]
    String DIR = SM + "." + ID

    File read_pbi = sub(reads_bam, ".bam$", ".bam.pbi")
    call ShardLongReads as t_04_ShardLongReads {
        input:
            unaligned_bam = reads_bam,
            unaligned_pbi = read_pbi,
            prefix = SM + "_shard",
            num_shards = primary_scatter_width,
    }

    call FindCCSReport as t_05_FindCCSReport {
        input:
            gcs_input_dir = gcs_input_dir
    }

    scatter (sharded_reads in t_04_ShardLongReads.unmapped_shards) {

        #######################################################################
        ## Setup
        #######################################################################

        String fbmrq_prefix = basename(sharded_reads, ".bam")

        # Filter out the kinetics tags from PB files:
        call RemoveKineticsTags as t_06_RemoveKineticsTags {
            input:
                bam = sharded_reads,
                prefix = SM + "_kinetics_removed"
        }

        #######################################################################
        ## Handle CCS Corrected Reads
        #######################################################################

        # 1 - filter the reads by the minimum read quality:
        call Bamtools as t_07_FilterByMinQual {
            input:
                bamfile = t_06_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_good_reads",
                cmd = "filter",
                args = '-tag "rq":">=' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # Annotate our CCS Corrected reads:
        call Annotate as t_08_LongbowAnnotateCCSReads {
            input:
                reads = t_07_FilterByMinQual.bam_out,
                model = mas_seq_model
        }

        # 6: Longbow filter ccs reclaimable reads
        call Filter as t_09_FilterCCSReads {
            input:
                bam = t_08_LongbowAnnotateCCSReads.annotated_bam,
                prefix = SM + "_ccs_corrected_subshard",
                model = mas_seq_model
        }

        call PBIndex as t_10_PbIndexLongbowAnnotatedCCSPassedReads {
            input:
                bam = t_09_FilterCCSReads.passed_reads
        }

        # Shard these reads even wider so we can make sure we don't run out of memory:
        call ShardLongReads as t_11_ShardCorrectedReads {
            input:
                unaligned_bam = t_09_FilterCCSReads.passed_reads,
                unaligned_pbi = t_10_PbIndexLongbowAnnotatedCCSPassedReads.pbindex,
                prefix = SM + "_ccs_corrected_longbow_annotated_subshard",
                num_shards = secondary_scatter_width,
        }

        # Segment our arrays into individual array elements:
        scatter (corrected_shard in t_11_ShardCorrectedReads.unmapped_shards) {
            call Segment as t_12_SegmentCCSAnnotatedReads {
                input:
                    annotated_reads = corrected_shard,
                    model = mas_seq_model
            }
        }

        # Merge all outputs of Longbow Annotate / Segment:
        call MergeBams as t_13_MergeCCSArrayElements_1 {
            input:
                bams = t_12_SegmentCCSAnnotatedReads.segmented_bam,
                prefix = SM + "_ccs_corrected_ArrayElements_intermediate_1"
        }

        # Align our array elements to the transcriptome:
        call Minimap2 as t_14_AlignCCSArrayElements {
            input:
                reads      = [ t_13_MergeCCSArrayElements_1.merged_bam ],
                ref_fasta  = transcriptome_ref_fasta,
                map_preset = "asm20"
        }

        # Align our array elements to the genome in splice-aware mode:
        call Minimap2 as t_15_AlignCCSArrayElementsToGenome {
            input:
                reads      = [ t_13_MergeCCSArrayElements_1.merged_bam ],
                ref_fasta  = ref_fasta,
                map_preset = "splice:hq"
        }

        # We need to restore the annotations we created with the 10x tool to the transcriptome aligned reads.
        call RestoreAnnotationstoAlignedBam as t_16_RestoreAnnotationsToTranscriptomeAlignedCCSBam {
            input:
                annotated_bam_file = t_13_MergeCCSArrayElements_1.merged_bam,
                aligned_bam_file = t_14_AlignCCSArrayElements.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # We need to restore the annotations we created with the 10x tool to the genome aligned reads.
        call RestoreAnnotationstoAlignedBam as t_17_RestoreAnnotationsToGenomeAlignedCCSBam {
            input:
                annotated_bam_file = t_13_MergeCCSArrayElements_1.merged_bam,
                aligned_bam_file = t_15_AlignCCSArrayElementsToGenome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # To properly count our transcripts we must throw away the non-primary and unaligned reads:
        # Filter out all non-primary transcriptome alignments:
        call FilterReadsBySamFlags as t_18_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads {
            input:
                bam = t_16_RestoreAnnotationsToTranscriptomeAlignedCCSBam.output_bam,
                sam_flags = "2308",
                prefix = SM + "_ccs_corrected_ArrayElements_Annotated_Transcriptome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }
        # Filter out all non-primary genome alignments:
        call FilterReadsBySamFlags as t_19_RemoveUnmappedAndNonPrimaryCCSGenomeReads {
            input:
                bam = t_17_RestoreAnnotationsToGenomeAlignedCCSBam.output_bam,
                sam_flags = "2308",
                prefix = SM + "_ccs_corrected_ArrayElements_Annotated_Genome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }

        #######################################################################
        ## Handle CCS Uncorrected Reads
        #######################################################################

        # 1.5 - Get the "rejected" reads:
        call Bamtools as t_20_GetCcsRejectedReads {
            input:
                bamfile = t_06_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_rejected_reads",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # 2 - Get reads we can reclaim:
        call Bamtools as t_21_ExtractCcsReclaimableReads {
            input:
                bamfile = t_06_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_reads_for_ccs_reclamation",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '" -length "<=' + max_reclamation_length + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # Annotate our CCS uncorrected (reclaimable) reads
        call Annotate as t_22_AnnotateReclaimableReads {
            input:
                reads = t_21_ExtractCcsReclaimableReads.bam_out,
                model = mas_seq_model
        }

        call PBIndex as t_23_PbIndexLongbowAnnotatedReclaimedReads {
            input:
                bam = t_22_AnnotateReclaimableReads.annotated_bam
        }

        # 6: Longbow filter ccs reclaimable reads
        call Filter as t_24_FilterReclaimableReads {
            input:
                bam = t_22_AnnotateReclaimableReads.annotated_bam,
                bam_pbi = t_23_PbIndexLongbowAnnotatedReclaimedReads.pbindex,
                prefix = SM + "_subshard",
                model = mas_seq_model
        }

        call PBIndex as t_25_PbIndexLongbowAnnotatedReclaimedPassedReads {
            input:
                bam = t_24_FilterReclaimableReads.passed_reads
        }

        # Shard these reads even wider so we can make sure we don't run out of memory:
        call ShardLongReads as t_26_ShardReclaimedReads {
            input:
                unaligned_bam = t_24_FilterReclaimableReads.passed_reads,
                unaligned_pbi = t_25_PbIndexLongbowAnnotatedReclaimedPassedReads.pbindex,
                prefix = SM + "_ccs_reclaimed_longbow_annotated_subshard",
                num_shards = secondary_scatter_width,
        }

        # Segment our arrays into individual array elements:
        scatter (corrected_shard in t_26_ShardReclaimedReads.unmapped_shards) {
            call Segment as t_27_SegmentReclaimedAnnotatedReads {
                input:
                    annotated_reads = corrected_shard,
                    model = mas_seq_model
            }
        }

        # Merge all outputs of Longbow Annotate / Segment:
        call MergeBams as t_28_MergeReclaimedArrayElements_1 {
            input:
                bams = t_27_SegmentReclaimedAnnotatedReads.segmented_bam,
                prefix = SM + "_ccs_reclaimed_ArrayElements_intermediate_1"
        }

        # Align our array elements to the transcriptome:
        call Minimap2 as t_29_AlignReclaimedArrayElements {
            input:
                reads      = [ t_28_MergeReclaimedArrayElements_1.merged_bam ],
                ref_fasta  = transcriptome_ref_fasta,
                map_preset = "asm20"
        }

        # Align our array elements to the genome in splice-aware mode:
        call Minimap2 as t_30_AlignReclaimedArrayElementsToGenome {
            input:
                reads      = [ t_28_MergeReclaimedArrayElements_1.merged_bam ],
                ref_fasta  = ref_fasta,
                map_preset = "splice:hq"
        }

        # We need to restore the annotations we created with the 10x tool to the transcriptome aligned reads.
        call RestoreAnnotationstoAlignedBam as t_31_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam {
            input:
                annotated_bam_file = t_28_MergeReclaimedArrayElements_1.merged_bam,
                aligned_bam_file = t_29_AlignReclaimedArrayElements.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # We need to restore the annotations we created with the 10x tool to the genome aligned reads.
        call RestoreAnnotationstoAlignedBam as t_32_RestoreAnnotationsToGenomeAlignedReclaimedBam {
            input:
                annotated_bam_file = t_28_MergeReclaimedArrayElements_1.merged_bam,
                aligned_bam_file = t_30_AlignReclaimedArrayElementsToGenome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # To properly count our transcripts we must throw away the non-primary and unaligned reads:
        # Filter out all non-primary transcriptome alignments:
        call FilterReadsBySamFlags as t_33_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads {
            input:
                bam = t_31_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam.output_bam,
                sam_flags = "2308",
                prefix = SM + "_ccs_reclaimed_ArrayElements_Annotated_Transcriptome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }
        # Filter out all non-primary genome alignments:
        call FilterReadsBySamFlags as t_34_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads {
            input:
                bam = t_32_RestoreAnnotationsToGenomeAlignedReclaimedBam.output_bam,
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
    call MergeBams as t_35_MergeCCSReads { input: bams = t_07_FilterByMinQual.bam_out, prefix = SM + "_ccs_reads" }
    call MergeBams as t_36_MergeAnnotatedCCSReads { input: bams = t_08_LongbowAnnotateCCSReads.annotated_bam, prefix = SM + "_ccs_reads_annotated" }

    # Merge all CCS bams together for this Subread BAM:
    call MergeBams as t_37_MergePassedCCSReads { input: bams = t_09_FilterCCSReads.passed_reads, prefix = SM + "_ccs_longbow_passed_reads" }
    call MergeBams as t_38_MergeFailedCCSReads { input: bams = t_09_FilterCCSReads.failed_reads, prefix = SM + "_ccs_longbow_failed_reads" }
    call MergeBams as t_39_MergeCCSArrayElements { input: bams = t_13_MergeCCSArrayElements_1.merged_bam, prefix = SM + "_ccs_array_elements" }
    call MergeBams as t_40_MergeTranscriptomeAlignedCCSArrayElements { input: bams = t_16_RestoreAnnotationsToTranscriptomeAlignedCCSBam.output_bam, prefix = SM + "_ccs_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call MergeBams as t_41_MergeGenomeAlignedCCSArrayElements { input: bams = t_17_RestoreAnnotationsToGenomeAlignedCCSBam.output_bam, prefix = SM + "_ccs_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call MergeBams as t_42_MergePrimaryTranscriptomeAlignedCCSArrayElements { input: bams = t_18_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads.output_bam, prefix = SM + "_ccs_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call MergeBams as t_43_MergePrimaryGenomeAlignedCCSArrayElements { input: bams = t_19_RemoveUnmappedAndNonPrimaryCCSGenomeReads.output_bam, prefix = SM + "_ccs_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    #################################################
    # CCS Reclaimed:
    call MergeBams as t_44_MergeCCSReclaimableReads { input: bams = t_21_ExtractCcsReclaimableReads.bam_out, prefix = SM + "_reclaimable_reads" }
    call MergeBams as t_45_MergeAnnotatedCCSReclaimableReads { input: bams = t_22_AnnotateReclaimableReads.annotated_bam, prefix = SM + "_reclaimable_reads_annotated" }

    # Merge all CCS bams together for this Subread BAM:
    call MergeBams as t_46_MergePassedCCSReclaimedReads { input: bams = t_24_FilterReclaimableReads.passed_reads, prefix = SM + "_reclaimable_longbow_passed_reads" }
    call MergeBams as t_47_MergeFailedCCSReclaimableReads { input: bams = t_24_FilterReclaimableReads.failed_reads, prefix = SM + "_reclaimed_longbow_failed_reads" }
    call MergeBams as t_48_MergeCCSReclaimedArrayElements { input: bams = t_28_MergeReclaimedArrayElements_1.merged_bam, prefix = SM + "_reclaimed_array_elements" }
    call MergeBams as t_49_MergeTranscriptomeAlignedCCSReclaimedArrayElements { input: bams = t_31_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam.output_bam, prefix = SM + "_reclaimed_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call MergeBams as t_50_MergeGenomeAlignedCCSReclaimedArrayElements { input: bams = t_32_RestoreAnnotationsToGenomeAlignedReclaimedBam.output_bam, prefix = SM + "_reclaimed_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call MergeBams as t_51_MergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements { input: bams = t_33_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads.output_bam, prefix = SM + "_reclaimed_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call MergeBams as t_52_MergePrimaryGenomeAlignedCCSReclaimedArrayElements { input: bams = t_34_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads.output_bam, prefix = SM + "_reclaimed_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    #################################################
    # Overall:
    call MergeBams as t_53_MergeAllAnnotatedReads { input: bams = flatten([t_08_LongbowAnnotateCCSReads.annotated_bam, t_22_AnnotateReclaimableReads.annotated_bam]), prefix = SM + "_all_annotated_reads" }
    call MergeBams as t_54_MergeAllLongbowPassedReads { input: bams = flatten([t_09_FilterCCSReads.passed_reads, t_24_FilterReclaimableReads.passed_reads]), prefix = SM + "_all_longbow_passed_reads" }
    call MergeBams as t_55_MergeAllLongbowFailedReads { input: bams = flatten([t_09_FilterCCSReads.failed_reads, t_24_FilterReclaimableReads.failed_reads]), prefix = SM + "_all_longbow_failed_reads" }

    call MergeBams as t_56_MergeAllArrayElements { input: bams = flatten([t_13_MergeCCSArrayElements_1.merged_bam, t_28_MergeReclaimedArrayElements_1.merged_bam]), prefix = SM + "_all_array_elements", runtime_attr_override = merge_extra_cpu_attrs }
    call MergeBams as t_57_MergeAllTranscriptomeAlignedArrayElements { input: bams = flatten([t_16_RestoreAnnotationsToTranscriptomeAlignedCCSBam.output_bam, t_31_RestoreAnnotationsToTranscriptomeAlignedReclaimedBam.output_bam]), prefix = SM + "_all_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call MergeBams as t_58_MergeAllGenomeAlignedArrayElements { input: bams = flatten([t_17_RestoreAnnotationsToGenomeAlignedCCSBam.output_bam, t_32_RestoreAnnotationsToGenomeAlignedReclaimedBam.output_bam]), prefix = SM + "_all_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call MergeBams as t_59_MergeAllPrimaryTranscriptomeAlignedArrayElements { input: bams = flatten([t_18_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads.output_bam, t_33_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads.output_bam]), prefix = SM + "_all_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call MergeBams as t_60_MergeAllPrimaryGenomeAlignedArrayElements { input: bams = flatten([t_19_RemoveUnmappedAndNonPrimaryCCSGenomeReads.output_bam, t_34_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads.output_bam]), prefix = SM + "_all_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    #################################################
    #   ___      ____
    #  / _ \    / ___|
    # | | | |  | |
    # | |_| |  | |___
    #  \__\_\   \____|
    #
    #################################################

    # Get stats on CCS reads:
    call Stats as t_61_CCS_longbow_stats {
        input:
            reads = t_36_MergeAnnotatedCCSReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_CCS_Corrected"
    }

    # Get stats on Reclaimable reads:
    call Stats as t_62_Reclaimable_longbow_stats {
        input:
            reads = t_45_MergeAnnotatedCCSReclaimableReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_CCS_Reclaimable"
    }

    # Get stats on Reclaimed reads:
    call Stats as t_63_Reclaimed_longbow_stats {
        input:
            reads = t_46_MergePassedCCSReclaimedReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_CCS_Reclaimed"
    }

    # Get stats on All reads (overall stats):
    call Stats as t_64_Passed_longbow_stats {
        input:
            reads = t_54_MergeAllLongbowPassedReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_All_Longbow_Passed"
    }

    # Get stats on All reads (overall stats):
    call Stats as t_65_Overall_longbow_stats {
        input:
            reads = t_53_MergeAllAnnotatedReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_Overall"
    }

    ######################################################################
    #     ___        _               _
    #    / _ \ _   _| |_ _ __  _   _| |_
    #   | | | | | | | __| '_ \| | | | __|
    #   | |_| | |_| | |_| |_) | |_| | |_
    #    \___/ \__,_|\__| .__/ \__,_|\__|
    #                   |_|
    ######################################################################

    output {
        File ccs_report = t_05_FindCCSReport.ccs_report

        File ccs_corrected_reads = t_35_MergeCCSReads.merged_bam
        File ccs_corrected_reads_index = t_35_MergeCCSReads.merged_bai
        File annotated_ccs_corrected_reads = t_36_MergeAnnotatedCCSReads.merged_bam
        File annotated_ccs_corrected_reads_index = t_36_MergeAnnotatedCCSReads.merged_bai

        File ccs_reclaimable_reads = t_44_MergeCCSReclaimableReads.merged_bam
        File ccs_reclaimable_reads_index = t_44_MergeCCSReclaimableReads.merged_bai
        File annotated_ccs_reclaimable_reads = t_45_MergeAnnotatedCCSReclaimableReads.merged_bam
        File annotated_ccs_reclaimable_reads_index = t_45_MergeAnnotatedCCSReclaimableReads.merged_bai
        File ccs_reclaimed_reads = t_46_MergePassedCCSReclaimedReads.merged_bam
        File ccs_reclaimed_reads_index = t_46_MergePassedCCSReclaimedReads.merged_bai
        File ccs_unreclaimable_reads = t_47_MergeFailedCCSReclaimableReads.merged_bam
        File ccs_unreclaimable_reads_index = t_47_MergeFailedCCSReclaimableReads.merged_bai

        File all_annotated_reads = t_53_MergeAllAnnotatedReads.merged_bam
        File all_annotated_reads_index = t_53_MergeAllAnnotatedReads.merged_bai
        File all_longbow_passed_reads = t_54_MergeAllLongbowPassedReads.merged_bam
        File all_longbow_passed_reads_index = t_54_MergeAllLongbowPassedReads.merged_bai
        File all_longbow_failed_reads = t_55_MergeAllLongbowFailedReads.merged_bam
        File all_longbow_failed_reads_index = t_55_MergeAllLongbowFailedReads.merged_bai

        File ccs_corrected_array_elements = t_39_MergeCCSArrayElements.merged_bam
        File ccs_corrected_array_elements_index = t_39_MergeCCSArrayElements.merged_bai
        File ccs_corrected_array_elements_tx_aligned = t_40_MergeTranscriptomeAlignedCCSArrayElements.merged_bam
        File ccs_corrected_array_elements_tx_aligned_index = t_40_MergeTranscriptomeAlignedCCSArrayElements.merged_bai
        File ccs_corrected_array_elements_genome_aligned = t_41_MergeGenomeAlignedCCSArrayElements.merged_bam
        File ccs_corrected_array_elements_genome_aligned_index = t_41_MergeGenomeAlignedCCSArrayElements.merged_bai
        File ccs_corrected_array_elements_tx_aligned_primary_alignments = t_42_MergePrimaryTranscriptomeAlignedCCSArrayElements.merged_bam
        File ccs_corrected_array_elements_tx_aligned_primary_alignments_index = t_42_MergePrimaryTranscriptomeAlignedCCSArrayElements.merged_bai
        File ccs_corrected_array_elements_genome_aligned_primary_alignments = t_43_MergePrimaryGenomeAlignedCCSArrayElements.merged_bam
        File ccs_corrected_array_elements_genome_aligned_primary_alignments_index = t_43_MergePrimaryGenomeAlignedCCSArrayElements.merged_bai

        File ccs_reclaimed_array_elements  = t_48_MergeCCSReclaimedArrayElements.merged_bam
        File ccs_reclaimed_array_elements_index  = t_48_MergeCCSReclaimedArrayElements.merged_bai
        File ccs_reclaimed_array_elements_tx_aligned  = t_49_MergeTranscriptomeAlignedCCSReclaimedArrayElements.merged_bam
        File ccs_reclaimed_array_elements_tx_aligned_index  = t_49_MergeTranscriptomeAlignedCCSReclaimedArrayElements.merged_bai
        File ccs_reclaimed_array_elements_genome_aligned  = t_50_MergeGenomeAlignedCCSReclaimedArrayElements.merged_bam
        File ccs_reclaimed_array_elements_genome_aligned_index  = t_50_MergeGenomeAlignedCCSReclaimedArrayElements.merged_bai
        File ccs_reclaimed_array_elements_tx_aligned_primary_alignments  = t_51_MergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements.merged_bam
        File ccs_reclaimed_array_elements_tx_aligned_primary_alignments_index  = t_51_MergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements.merged_bai
        File ccs_reclaimed_array_elements_genome_aligned_primary_alignments  = t_52_MergePrimaryGenomeAlignedCCSReclaimedArrayElements.merged_bam
        File ccs_reclaimed_array_elements_genome_aligned_primary_alignments_index  = t_52_MergePrimaryGenomeAlignedCCSReclaimedArrayElements.merged_bai

        File all_array_elements = t_56_MergeAllArrayElements.merged_bam
        File all_array_elements_index = t_56_MergeAllArrayElements.merged_bai
        File all_longbow_passed_array_elements_tx_aligned = t_57_MergeAllTranscriptomeAlignedArrayElements.merged_bam
        File all_longbow_passed_array_elements_tx_aligned_index = t_57_MergeAllTranscriptomeAlignedArrayElements.merged_bai
        File all_longbow_passed_array_elements_genome_aligned = t_58_MergeAllGenomeAlignedArrayElements.merged_bam
        File all_longbow_passed_array_elements_genome_aligned_index = t_58_MergeAllGenomeAlignedArrayElements.merged_bai
        File all_longbow_passed_array_elements_tx_aligned_primary_alignments = t_59_MergeAllPrimaryTranscriptomeAlignedArrayElements.merged_bam
        File all_longbow_passed_array_elements_tx_aligned_primary_alignments_index = t_59_MergeAllPrimaryTranscriptomeAlignedArrayElements.merged_bai
        File all_longbow_passed_array_elements_genome_aligned_primary_alignments = t_60_MergeAllPrimaryGenomeAlignedArrayElements.merged_bam
        File all_longbow_passed_array_elements_genome_aligned_primary_alignments_index = t_60_MergeAllPrimaryGenomeAlignedArrayElements.merged_bai

        File ccs_corrected_longbow_stats_summary_stats = t_61_CCS_longbow_stats.summary_stats
        File ccs_corrected_longbow_stats_array_length_counts_plot_png = t_61_CCS_longbow_stats.array_length_counts_plot_png
        File ccs_corrected_longbow_stats_array_length_counts_plot_svg = t_61_CCS_longbow_stats.array_length_counts_plot_svg
        File ccs_corrected_longbow_stats_ligation_heatmap_nn_png = t_61_CCS_longbow_stats.ligation_heatmap_nn_png
        File ccs_corrected_longbow_stats_ligation_heatmap_nn_svg = t_61_CCS_longbow_stats.ligation_heatmap_nn_svg
        File ccs_corrected_longbow_stats_ligation_heatmap_png = t_61_CCS_longbow_stats.ligation_heatmap_png
        File ccs_corrected_longbow_stats_ligation_heatmap_svg = t_61_CCS_longbow_stats.ligation_heatmap_svg
        File ccs_corrected_longbow_stats_ligation_heatmap_nn_reduced_png = t_61_CCS_longbow_stats.ligation_heatmap_nn_reduced_png
        File ccs_corrected_longbow_stats_ligation_heatmap_nn_reduced_svg = t_61_CCS_longbow_stats.ligation_heatmap_nn_reduced_svg
        File ccs_corrected_longbow_stats_ligation_heatmap_reduced_png = t_61_CCS_longbow_stats.ligation_heatmap_reduced_png
        File ccs_corrected_longbow_stats_ligation_heatmap_reduced_svg = t_61_CCS_longbow_stats.ligation_heatmap_reduced_svg

        File ccs_reclaimable_longbow_stats_summary_stats = t_62_Reclaimable_longbow_stats.summary_stats
        File ccs_reclaimable_longbow_stats_array_length_counts_plot_png = t_62_Reclaimable_longbow_stats.array_length_counts_plot_png
        File ccs_reclaimable_longbow_stats_array_length_counts_plot_svg = t_62_Reclaimable_longbow_stats.array_length_counts_plot_svg
        File ccs_reclaimable_longbow_stats_ligation_heatmap_nn_png = t_62_Reclaimable_longbow_stats.ligation_heatmap_nn_png
        File ccs_reclaimable_longbow_stats_ligation_heatmap_nn_svg = t_62_Reclaimable_longbow_stats.ligation_heatmap_nn_svg
        File ccs_reclaimable_longbow_stats_ligation_heatmap_png = t_62_Reclaimable_longbow_stats.ligation_heatmap_png
        File ccs_reclaimable_longbow_stats_ligation_heatmap_svg = t_62_Reclaimable_longbow_stats.ligation_heatmap_svg
        File ccs_reclaimable_longbow_stats_ligation_heatmap_nn_reduced_png = t_62_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_png
        File ccs_reclaimable_longbow_stats_ligation_heatmap_nn_reduced_svg = t_62_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_svg
        File ccs_reclaimable_longbow_stats_ligation_heatmap_reduced_png = t_62_Reclaimable_longbow_stats.ligation_heatmap_reduced_png
        File ccs_reclaimable_longbow_stats_ligation_heatmap_reduced_svg = t_62_Reclaimable_longbow_stats.ligation_heatmap_reduced_svg

        File ccs_reclaimed_longbow_stats_summary_stats = t_63_Reclaimed_longbow_stats.summary_stats
        File ccs_reclaimed_longbow_stats_array_length_counts_plot_png = t_63_Reclaimed_longbow_stats.array_length_counts_plot_png
        File ccs_reclaimed_longbow_stats_array_length_counts_plot_svg = t_63_Reclaimed_longbow_stats.array_length_counts_plot_svg
        File ccs_reclaimed_longbow_stats_ligation_heatmap_nn_png = t_63_Reclaimed_longbow_stats.ligation_heatmap_nn_png
        File ccs_reclaimed_longbow_stats_ligation_heatmap_nn_svg = t_63_Reclaimed_longbow_stats.ligation_heatmap_nn_svg
        File ccs_reclaimed_longbow_stats_ligation_heatmap_png = t_63_Reclaimed_longbow_stats.ligation_heatmap_png
        File ccs_reclaimed_longbow_stats_ligation_heatmap_svg = t_63_Reclaimed_longbow_stats.ligation_heatmap_svg
        File ccs_reclaimed_longbow_stats_ligation_heatmap_nn_reduced_png = t_63_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_png
        File ccs_reclaimed_longbow_stats_ligation_heatmap_nn_reduced_svg = t_63_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_svg
        File ccs_reclaimed_longbow_stats_ligation_heatmap_reduced_png = t_63_Reclaimed_longbow_stats.ligation_heatmap_reduced_png
        File ccs_reclaimed_longbow_stats_ligation_heatmap_reduced_svg = t_63_Reclaimed_longbow_stats.ligation_heatmap_reduced_svg

        File overall_longbow_stats_summary_stats = t_65_Overall_longbow_stats.summary_stats
        File overall_longbow_stats_array_length_counts_plot_png = t_65_Overall_longbow_stats.array_length_counts_plot_png
        File overall_longbow_stats_array_length_counts_plot_svg = t_65_Overall_longbow_stats.array_length_counts_plot_svg
        File overall_longbow_stats_ligation_heatmap_nn_png = t_65_Overall_longbow_stats.ligation_heatmap_nn_png
        File overall_longbow_stats_ligation_heatmap_nn_svg = t_65_Overall_longbow_stats.ligation_heatmap_nn_svg
        File overall_longbow_stats_ligation_heatmap_png = t_65_Overall_longbow_stats.ligation_heatmap_png
        File overall_longbow_stats_ligation_heatmap_svg = t_65_Overall_longbow_stats.ligation_heatmap_svg
        File overall_longbow_stats_ligation_heatmap_nn_reduced_png = t_65_Overall_longbow_stats.ligation_heatmap_nn_reduced_png
        File overall_longbow_stats_ligation_heatmap_nn_reduced_svg = t_65_Overall_longbow_stats.ligation_heatmap_nn_reduced_svg
        File overall_longbow_stats_ligation_heatmap_reduced_png = t_65_Overall_longbow_stats.ligation_heatmap_reduced_png
        File overall_longbow_stats_ligation_heatmap_reduced_svg = t_65_Overall_longbow_stats.ligation_heatmap_reduced_svg

        File longbow_passed_longbow_stats_summary_stats = t_64_Passed_longbow_stats.summary_stats
        File longbow_passed_longbow_stats_array_length_counts_plot_png = t_64_Passed_longbow_stats.array_length_counts_plot_png
        File longbow_passed_longbow_stats_array_length_counts_plot_svg = t_64_Passed_longbow_stats.array_length_counts_plot_svg
        File longbow_passed_longbow_stats_ligation_heatmap_nn_png = t_64_Passed_longbow_stats.ligation_heatmap_nn_png
        File longbow_passed_longbow_stats_ligation_heatmap_nn_svg = t_64_Passed_longbow_stats.ligation_heatmap_nn_svg
        File longbow_passed_longbow_stats_ligation_heatmap_png = t_64_Passed_longbow_stats.ligation_heatmap_png
        File longbow_passed_longbow_stats_ligation_heatmap_svg = t_64_Passed_longbow_stats.ligation_heatmap_svg
        File longbow_passed_longbow_stats_ligation_heatmap_nn_reduced_png = t_64_Passed_longbow_stats.ligation_heatmap_nn_reduced_png
        File longbow_passed_longbow_stats_ligation_heatmap_nn_reduced_svg = t_64_Passed_longbow_stats.ligation_heatmap_nn_reduced_svg
        File longbow_passed_longbow_stats_ligation_heatmap_reduced_png = t_64_Passed_longbow_stats.ligation_heatmap_reduced_png
        File longbow_passed_longbow_stats_ligation_heatmap_reduced_svg = t_64_Passed_longbow_stats.ligation_heatmap_reduced_svg
    }
}

########################################################################################################################
########################################################################################################################
########################################################################################################################
#        _____         _
#       |_   _|_ _ ___| | _____
#         | |/ _` / __| |/ / __|
#         | | (_| \__ \   <\__ \
#         |_|\__,_|___/_|\_\___/
#
########################################################################################################################
########################################################################################################################
########################################################################################################################

    ######################################################################
    #     ____  ____  _   _ _   _ _
    #    |  _ \| __ )| | | | |_(_) |___
    #    | |_) |  _ \| | | | __| | / __|
    #    |  __/| |_) | |_| | |_| | \__ \
    #    |_|   |____/ \___/ \__|_|_|___/
    #
    ######################################################################

task FindBams {
    input {
        String gcs_input_dir

        RuntimeAttr? runtime_attr_override
    }

    String indir = sub(gcs_input_dir, "/$", "")

    command <<<
        # Commenting out this line to make it able to run even if there are no subreads / reads files:
        # set -euxo pipefail

        # Get our file lists here:
        gsutil ls ~{indir}/**subreads.bam | sort > subread_bams.txt
        gsutil ls ~{indir}/**.reads.bam | sort > ccs_bams.txt

        # Check our file lists to see which are populated:
        if [[ $( wc -l subread_bams.txt | awk '{print $1}' ) -gt 0 ]] ; then
            echo 'Found subreads!'
            echo "true" > has_subreads.txt
        else
            echo 'Did not find subreads!'
            echo "false" > has_subreads.txt
        fi

        if [[ $( wc -l ccs_bams.txt | awk '{print $1}' ) -gt 0 ]] ; then
            echo 'Found CCS reads!'
            echo "true" > has_ccs_reads.txt
        else
            echo 'Did not find CCS reads!'
            echo "false" > has_ccs_reads.txt
        fi
    >>>

    output {
        Array[String] subread_bams = read_lines("subread_bams.txt")
        Array[String] ccs_bams = read_lines("ccs_bams.txt")

        Boolean has_subreads = read_boolean("has_subreads.txt")
        Boolean has_ccs_reads = read_boolean("has_ccs_reads.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
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

task FindCCSReport {
    input {
        String gcs_input_dir

        RuntimeAttr? runtime_attr_override
    }

    String indir = sub(gcs_input_dir, "/$", "")

    command <<<
        set -euxo pipefail
        gsutil ls ~{indir}/**.ccs_reports.txt | sort > ccs_reports.txt
    >>>

    output {
        # There should be only one CCS report here:
        File ccs_report = read_lines("ccs_reports.txt")[0]
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
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

task GetRunInfo {
    input {
        String subread_bam

        String? bam_suffix

        RuntimeAttr? runtime_attr_override
    }

    String gcs_dir = sub(subread_bam, basename(subread_bam), "")

    String bam_suffix_arg = if defined(bam_suffix) then " --BS " else ""

    command <<<
        set -x

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

        # We need to update detect_run_info.py to make it sanitize fields.
        # The `sed` statement here is a hack to get around an issue.
        python /usr/local/bin/detect_run_info.py ~{gcs_dir} ~{bam_suffix_arg}~{default="" sep=" --BS " bam_suffix} > run_info.txt
    >>>

    output {
        Map[String, String] run_info = read_map("run_info.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.30"
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

task RemoveKineticsTags {
    input {
        File bam
        String prefix = "kinetics_removed"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(bam, "GB"))
    command <<<
        set -x

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        samtools view -hb -@$np -x fi -x ri -x fp -x rp -x fn -x rn ~{bam} > ~{prefix}.bam
    >>>

    output {
        File bam_file = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.21"
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

task ShardLongReads {
    input {
        File unaligned_bam
        File unaligned_pbi

        Int num_shards = 300
        Int num_threads = 8

        String prefix = "shard"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3*ceil(size(unaligned_bam, "GB") + size(unaligned_pbi, "GB"))
    Int mem = ceil(25*size(unaligned_pbi, "MB")/1000)

    command <<<
        set -x

        python3 /usr/local/bin/shard_bam.py \
            -n ~{num_shards} \
            -t ~{num_threads} \
            -i ~{unaligned_pbi} \
            -p ~{prefix} \
            ~{unaligned_bam}
    >>>

    output {
        Array[File] unmapped_shards = glob("~{prefix}*.bam")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_threads,
        mem_gb:             mem,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.21"
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

task PBIndex {
    input {
        File bam
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB"))

    String base_name = basename(bam)

    command <<<
        set -euxo pipefail

        # Run PBIndex:
        pbindex ~{bam}
        mv ~{bam}.pbi ~{base_name}.pbi
    >>>

    output {
        File pbindex = "~{base_name}.pbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,             # 1 preemptible try here so that we have even money with our local SSDs below.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.21"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"  # LOCAL here is a local SSD - much faster and even money with normal disk if preemptible
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

    ######################################################################
    #     _   _ _   _ _
    #    | | | | |_(_) |___
    #    | | | | __| | / __|
    #    | |_| | |_| | \__ \
    #     \___/ \__|_|_|___/
    #
    ######################################################################

task MergeBams {
    input {
        Array[File] bams
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bams:   "input array of BAMs to be merged"
        prefix: "[default-valued] prefix for output BAM"
    }

    Int disk_size = 4*ceil(size(bams, "GB"))

    command <<<
        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        samtools merge -p -c -@$np --no-PG ~{prefix}.bam ~{sep=" " bams}
        samtools index ~{prefix}.bam
    >>>

    output {
        File merged_bam = "~{prefix}.bam"
        File merged_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             20,
        disk_gb:            disk_size,
        disk_type:          "HDD",
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
task Bamtools {
    input {
        File bamfile
        String cmd
        String args

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 1 + ceil(2 * size(bamfile, "GiB"))

    command <<<
        bamtools ~{cmd} -in ~{bamfile} -out ~{prefix}.bam ~{args}
    >>>

    output {
        File bam_out = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             2,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9.beta"
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

task FailWithWarning {
    input {
        String warning
    }
    command <<<
        set -e

        echo "~{warning}"
        echo "~{warning}" 1>&2
        false
    >>>
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
    runtime {
        cpu:                    default_attr.cpu_cores
        memory:                 default_attr.mem_gb + " GiB"
        disks: "local-disk " +  default_attr.disk_gb + " HDD"
        bootDiskSizeGb:         default_attr.boot_disk_gb
        preemptible:            default_attr.preemptible_tries
        maxRetries:             default_attr.max_retries
        docker:                 default_attr.docker
    }
}
task FilterReadsBySamFlags {
    meta {
        description : "Filter reads based on sam flags.  Reads with ANY of the given flags will be removed from the given dataset."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File bam
        String sam_flags

        String prefix = "filtered_reads"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:   "BAM file to be filtered."
        sam_flags: "Flags for which to remove reads.  Reads with ANY of the given flags will be removed from the given dataset."
        prefix : "[Optional] Prefix string to name the output file (Default: filtered_reads)."
    }

    Int disk_size = 20 + ceil(11 * size(bam, "GiB"))

    command <<<

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        samtools view -h -b -F ~{sam_flags} -@$np ~{bam} > ~{prefix}.bam
        samtools index -@$np ~{prefix}.bam
    >>>

    output {
        File output_bam = "~{prefix}.bam"
        File output_bam_index = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
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

    ######################################################################
    #      _    _ _             ____                _
    #     / \  | (_) __ _ _ __ |  _ \ ___  __ _  __| |___
    #    / _ \ | | |/ _` | '_ \| |_) / _ \/ _` |/ _` / __|
    #   / ___ \| | | (_| | | | |  _ <  __/ (_| | (_| \__ \
    #  /_/   \_\_|_|\__, |_| |_|_| \_\___|\__,_|\__,_|___/
    #               |___/
    ######################################################################


task Minimap2 {
    input {
        Array[File] reads
        File ref_fasta

        String map_preset

        String RG = ""

        String prefix = "out"
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:      "query sequences to be mapped and aligned"
        ref_fasta:  "reference fasta"
        RG:         "[optional] read group information to be supplied to parameter '-R' (note that tabs should be input as '\t')"
        map_preset: "preset to be used for minimap2 parameter '-x'"
        prefix:     "[default-valued] prefix for output BAM"
    }

    # 10x for the decompressed file size
    # 2x for potential for FASTQ and SAM files (from file conversion).
    # 2x for extra "just in case" space.
    # +1 to handle small files
    Int disk_size = 1 + 10*2*2*ceil(size(reads, "GB") + size(ref_fasta, "GB"))

    Int cpus = 4
    Int mem = 30

    # This is a hack to fix the WDL parsing of ${} variables:
    String DOLLAR = "$"
    command <<<
        set -euxo pipefail

        rg_len=$(echo -n '~{RG}' | wc -c | awk '{print $NF}')
        if [[ $rg_len -ne 0 ]] ; then
            # Sometimes we have to sanitize our read groups:
            sanitized_read_group=$( echo "~{RG}" | sed -e 's# .*##g' | sed 's#\t.*##g' )

            echo "Original Read Group: ~{RG}"
            echo "Sanitized Read Group: $sanitized_read_group"

            MAP_PARAMS="-ayYL --MD --eqx -x ~{map_preset} -R $sanitized_read_group -t ~{cpus} ~{ref_fasta}"
        else
            MAP_PARAMS="-ayYL --MD --eqx -x ~{map_preset} -t ~{cpus} ~{ref_fasta}"
        fi
        FILE="~{reads[0]}"
        FILES="~{sep=' ' reads}"

        # We write to a SAM file before sorting and indexing because rarely, doing everything
        # in a single one-liner leads to a truncated file error and segmentation fault of unknown
        # origin.  Separating these commands requires more resources, but is more reliable overall.

        if [[ "$FILE" =~ \.fastq$ ]] || [[ "$FILE" =~ \.fq$ ]]; then
            cat $FILES | minimap2 $MAP_PARAMS - > tmp.sam
        elif [[ "$FILE" =~ \.fastq.gz$ ]] || [[ "$FILE" =~ \.fq.gz$ ]]; then
            zcat $FILES | minimap2 $MAP_PARAMS - > tmp.sam
        elif [[ "$FILE" =~ \.fasta$ ]] || [[ "$FILE" =~ \.fa$ ]]; then
            cat $FILES | python3 /usr/local/bin/cat_as_fastq.py | minimap2 $MAP_PARAMS - > tmp.sam
        elif [[ "$FILE" =~ \.fasta.gz$ ]] || [[ "$FILE" =~ \.fa.gz$ ]]; then
            zcat $FILES | python3 /usr/local/bin/cat_as_fastq.py | minimap2 $MAP_PARAMS - > tmp.sam
        elif [[ "$FILE" =~ \.bam$ ]]; then

            # samtools fastq takes only 1 file at a time so we need to merge them together:
            for f in "~{sep=' ' reads}" ; do
                samtools fastq "$f"
            done > tmp.fastq

            echo "Memory info:"
            cat /proc/meminfo
            echo ""

            minimap2 $MAP_PARAMS tmp.fastq > tmp.sam
        else
            echo "Did not understand file format for '$FILE'"
            exit 1
        fi

        samtools sort -@~{cpus} -m~{mem}G --no-PG -o ~{prefix}.bam tmp.sam
        samtools index ~{prefix}.bam
    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             mem,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
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


    ######################################################################
    #    _____             __  __  _____           _
    #   |_   _|__ _ __     \ \/ / |_   _|__   ___ | |
    #     | |/ _ \ '_ \     \  /    | |/ _ \ / _ \| |
    #     | |  __/ | | |    /  \    | | (_) | (_) | |
    #     |_|\___|_| |_|___/_/\_\___|_|\___/ \___/|_|
    #                 |_____|  |_____|
    ######################################################################


task RestoreAnnotationstoAlignedBam {
    input {
        File annotated_bam_file
        File aligned_bam_file

        Array[String] tags_to_ignore = ["RG"]

        Int? boot_disk_size_gb
        Int? cpu
        Int? disk_space_gb
        Int? mem_gb
        Int? preemptible_attempts
    }

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 32 * 1024

    Float reads_size_gb = size(annotated_bam_file, "GiB") + size(aligned_bam_file, "GiB")
    Int default_disk_space_gb = 10 * ceil((reads_size_gb * 10) + 20)

    Int default_boot_disk_size_gb = 100

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb

    String timing_output_file = "timingInformation.txt"

    String memory_log_file = "memory_use.txt"
    String output_name = basename(aligned_bam_file, ".bam") + ".AnnotationsRestored.bam"

    String ignore_tags_arg = if (length(tags_to_ignore) != 0 ) then "--ignore-tags " else ""

    command {

        # Set up memory logging daemon:
        MEM_LOG_INTERVAL_s=5
        DO_MEMORY_LOG=true
        while $DO_MEMORY_LOG ; do
            date
            date +%s
            cat /proc/meminfo
            sleep $MEM_LOG_INTERVAL_s
        done >> ~{memory_log_file} &
        mem_pid=$!

        set -e
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        source activate 10x_tool

        python3 /lrma/restore_annotations_to_aligned_bam.py \
            --bam ~{annotated_bam_file} \
            --aligned-bam ~{aligned_bam_file} \
            ~{ignore_tags_arg} ~{default="" sep=" " tags_to_ignore} \
            --out-name ~{output_name}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ~{timing_output_file}

        # Stop the memory daemon softly.  Then stop it hard if it's not cooperating:
        set +e
        DO_MEMORY_LOG=false
        sleep $(($MEM_LOG_INTERVAL_s  * 2))
        kill -0 $mem_pid &> /dev/null
        if [ $? -ne 0 ] ; then
            kill -9 $mem_pid
        fi

        # Get and compute timing information:
        set +e
        elapsedTime=""
        which bc &> /dev/null ; bcr=$?
        which python3 &> /dev/null ; python3r=$?
        which python &> /dev/null ; pythonr=$?
        if [[ $bcr -eq 0 ]] ; then elapsedTime=`echo "scale=6;$endTime - $startTime" | bc`;
        elif [[ $python3r -eq 0 ]] ; then elapsedTime=`python3 -c "print( $endTime - $startTime )"`;
        elif [[ $pythonr -eq 0 ]] ; then elapsedTime=`python -c "print( $endTime - $startTime )"`;
        fi
        echo "Elapsed Time: $elapsedTime" >> ~{timing_output_file}
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-10x:0.1.15"
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: 0
        cpu: select_first([cpu, 1])
    }
    output {
      File output_bam        = "~{output_name}"
      File timing_info       = "~{timing_output_file}"
      File memory_info       = "~{memory_log_file}"
    }
}


    ######################################################################
    #    _                      _
    #   | |    ___  _ __   __ _| |__   _____      __
    #   | |   / _ \| '_ \ / _` | '_ \ / _ \ \ /\ / /
    #   | |__| (_) | | | | (_| | |_) | (_) \ V  V /
    #   |_____\___/|_| |_|\__, |_.__/ \___/ \_/\_/
    #                     |___/
    ######################################################################

task Annotate
{
    input {
        File reads
        String model = "mas15"
        String prefix = "longbow_annotated"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow annotate -t8 --model ~{model} -v INFO ~{reads} -o ~{prefix}.bam
    >>>

    output {
        File annotated_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.4.1"
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

task Segment
{
    input {
        File annotated_reads
        String model = "mas15"
        String prefix = "longbow_segmented"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 15*ceil(size(annotated_reads, "GB")) + 20

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow segment --model ~{model} -v INFO -s ~{annotated_reads} -o ~{prefix}.bam
    >>>

    output {
        File segmented_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.4.1"
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

task Filter {
    input {
        File bam
        String model = "mas15"

        String prefix = "reads"


        File? bam_pbi

        RuntimeAttr? runtime_attr_override
    }

    String pbi_arg = if defined(bam_pbi) then " --pbi " else ""
    Int disk_size = ceil(4*ceil(size(bam, "GB")) + size(bam_pbi, "GB"))

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow filter \
            -v INFO \
            ~{pbi_arg}~{default="" sep=" --pbi " bam_pbi} \
            --model ~{model} \
            ~{bam} \
            -o ~{prefix} 2> >(tee longbow_filter_log.txt >&2) # Get log data from stderr and reprint to stderr
    >>>

    output {
        File passed_reads = "~{prefix}_longbow_filter_passed.bam"
        File failed_reads = "~{prefix}_longbow_filter_failed.bam"
        File log = "longbow_filter_log.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.4.1"
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

task Stats
{
    input {
        File reads
        String model = "mas15"
        String prefix = "longbow_stats"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow stats -s --model ~{model} -v INFO -o ~{prefix} ~{reads}
    >>>

    output {
        File summary_stats = "~{prefix}_summary_stats.txt"

        File array_length_counts_plot_png = "~{prefix}_00_MAS-seq_Array_Length_Counts_~{model}.png"
        File array_length_counts_plot_svg = "~{prefix}_00_MAS-seq_Array_Length_Counts_~{model}.svg"

        File ligation_heatmap_nn_png = "~{prefix}_01_MAS-seq_Ligations_~{model}_no_numbers.png"
        File ligation_heatmap_nn_svg = "~{prefix}_01_MAS-seq_Ligations_~{model}_no_numbers.svg"

        File ligation_heatmap_png = "~{prefix}_02_MAS-seq_Ligations_~{model}.png"
        File ligation_heatmap_svg = "~{prefix}_02_MAS-seq_Ligations_~{model}.svg"

        File ligation_heatmap_nn_reduced_png = "~{prefix}_03_MAS-seq_Ligations_~{model}_reduced_no_numbers.png"
        File ligation_heatmap_nn_reduced_svg = "~{prefix}_03_MAS-seq_Ligations_~{model}_reduced_no_numbers.svg"

        File ligation_heatmap_reduced_png = "~{prefix}_04_MAS-seq_Ligations_~{model}_reduced.png"
        File ligation_heatmap_reduced_svg = "~{prefix}_04_MAS-seq_Ligations_~{model}_reduced.svg"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.4.1"
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
