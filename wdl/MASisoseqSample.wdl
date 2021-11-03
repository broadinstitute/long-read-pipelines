version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF
import "tasks/AlignReads.wdl" as AR
import "tasks/Ten_X_Tool.wdl" as TENX
import "tasks/Longbow.wdl" as LONGBOW
import "tasks/Structs.wdl"

workflow MASseqSample {
    meta {
        description : "This workflow processes MAS-seq data.  It annotates, segments, aligns, and (optionally) processes aligned reads into count matrices for downstream gene expression analysis.  For more information on MAS-seq see the preprint on biorXiv: https://www.biorxiv.org/content/10.1101/2021.10.01.462818"
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        Array[File] ccs_bams
        Array[File] ccs_pbis

        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/"

        # NOTE: Reference for un-split CCS reads:
        File genome_ref_map_file
        File transcriptome_ref_map_file

        File cell_barcode_whitelist

        File head_adapter_fasta
        File tail_adapter_fasta

        String sample_name

        String mas_seq_model = "mas15"

        Boolean is_SIRV_data = false
        Boolean generate_count_matrix = false

        Boolean use_model_for_barcode_annotation = false

        # Default here is 0 because ccs uncorrected reads all seem to have RQ = -1.
        # All pathologically long reads also have RQ = -1.
        # This way we preserve the vast majority of the data, even if it has low quality.
        # We can filter it out at later steps.
        Float min_read_quality = 0.0
        Int max_reclamation_length = 60000

        # Set up some meta parameters here so we can adjust for when we want things to go VERY fast:
        Int primary_scatter_width = 50
        Int secondary_scatter_width = 10

        String longbow_docker_version = "us.gcr.io/broad-dsp-lrma/lr-longbow:0.4.3"
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

    # Setup our reference information:
    Map[String, String] genome_ref_map = read_map(genome_ref_map_file)
    Map[String, String] transcriptome_ref_map = read_map(transcriptome_ref_map_file)

    # Make sure we have a good output directory:
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/MASseqSample/~{sample_name}/" + t_01_WdlExecutionStartTimestamp.timestamp_string

    # Set up our runtime attributes:
    RuntimeAttr disable_preemption_runtime_attrs = object {
        preemptible_tries: 0
    }

    RuntimeAttr fast_sharding_runtime_attrs = object {
        disk_type: "LOCAL",
        preemptible_tries: 0
    }

    RuntimeAttr filterReadsAttrs = object {
        cpu_cores: 4,
        preemptible_tries: 0
    }

    RuntimeAttr new_longbow_attrs = object {
        docker: longbow_docker_version
    }

    RuntimeAttr merge_extra_cpu_attrs = object {
        cpu_cores: 4
    }

    # Do some error checking before we begin:
    Int num_ccs_bams = length(ccs_bams)
    Int num_ccs_pbis = length(ccs_pbis)
    if (num_ccs_bams != num_ccs_pbis) {
        call Utils.FailWithWarning as t_02_WARN1 { input: warning = "Error: Num CCS bams is not the same as the number of CCS pbindex files  (~{num_ccs_bams} != ~{num_ccs_pbis}).  Aborting." }
    }

    String segment_extra_args = if use_model_for_barcode_annotation then "-b" else ""

    # Loop through all of our base input files and process them individually.
    # We'll merge everything together at the end.
    scatter ( i in range(length(ccs_bams))) {
        File ccs_bam = ccs_bams[i]
        File ccs_pbi = ccs_pbis[i]

        call PB.ShardLongReads as t_03_ShardLongReads {
            input:
                unaligned_bam = ccs_bam,
                unaligned_pbi = ccs_pbi,
                prefix = sample_name + "_shard",
                num_shards = 300,
                runtime_attr_override = fast_sharding_runtime_attrs,
        }

        scatter (sharded_reads in t_03_ShardLongReads.unmapped_shards) {

            #######################################################################
            ## Setup
            #######################################################################

            String fbmrq_prefix = basename(ccs_bam, ".bam")

            # Filter out the kinetics tags from PB files:
            call PB.RemoveKineticsTags as t_04_RemoveKineticsTags {
                input:
                    bam = ccs_bam,
                    prefix = sample_name + "_kinetics_removed"
            }

            #######################################################################
            ## Segment CCS Corrected Reads
            #######################################################################

            # 1 - filter the reads by the minimum read quality:
            call Utils.Bamtools as t_05_FilterByMinQual {
                input:
                    bamfile = t_04_RemoveKineticsTags.bam_file,
                    prefix = fbmrq_prefix + "_good_reads",
                    cmd = "filter",
                    args = '-tag "rq":">=' + min_read_quality + '"',
                    runtime_attr_override = disable_preemption_runtime_attrs
            }

            # Annotate our CCS Corrected reads:
            call LONGBOW.Annotate as t_06_LongbowAnnotateCCSReads {
                input:
                    reads = t_05_FilterByMinQual.bam_out,
                    model = mas_seq_model,
                    runtime_attr_override = new_longbow_attrs,
            }

            # 6: Longbow filter ccs reclaimable reads
            call LONGBOW.Filter as t_07_FilterCCSReads {
                input:
                    bam = t_06_LongbowAnnotateCCSReads.annotated_bam,
                    prefix = sample_name + "_ccs_corrected_subshard",
                    model = mas_seq_model,
                    runtime_attr_override = new_longbow_attrs,
            }

            call PB.PBIndex as t_08_PbIndexLongbowAnnotatedCCSPassedReads {
                input:
                    bam = t_07_FilterCCSReads.passed_reads
            }

            # Shard these reads even wider so we can make sure we don't run out of memory:
            call PB.ShardLongReads as t_09_ShardCorrectedReads {
                input:
                    unaligned_bam = t_07_FilterCCSReads.passed_reads,
                    unaligned_pbi = t_08_PbIndexLongbowAnnotatedCCSPassedReads.pbi,
                    prefix = sample_name + "_ccs_corrected_longbow_annotated_subshard",
                    num_shards = secondary_scatter_width,
            }

            # Segment our arrays into individual array elements:
            scatter (corrected_shard in t_09_ShardCorrectedReads.unmapped_shards) {
                call LONGBOW.Segment as t_10_SegmentCCSAnnotatedReads {
                    input:
                        annotated_reads = corrected_shard,
                        model = mas_seq_model,
                        runtime_attr_override = new_longbow_attrs,
                        extra_args = segment_extra_args
                }
            }

            # Merge all outputs of Longbow Annotate / Segment:
            call Utils.MergeBams as t_11_MergeCCSArrayElements_1 {
                input:
                    bams = t_10_SegmentCCSAnnotatedReads.segmented_bam,
                    prefix = sample_name + "_ccs_corrected_ArrayElements_intermediate_1"
            }

            #######################################################################
            ## Segment CCS Uncorrected Reads
            #######################################################################

            # 1.5 - Get the "rejected" reads:
            call Utils.Bamtools as t_16_GetCcsRejectedReads {
                input:
                    bamfile = t_04_RemoveKineticsTags.bam_file,
                    prefix = fbmrq_prefix + "_rejected_reads",
                    cmd = "filter",
                    args = '-tag "rq":"<' + min_read_quality + '"',
                    runtime_attr_override = disable_preemption_runtime_attrs
            }

            # 2 - Get reads we can reclaim:
            call Utils.Bamtools as t_17_ExtractCcsReclaimableReads {
                input:
                    bamfile = t_04_RemoveKineticsTags.bam_file,
                    prefix = fbmrq_prefix + "_reads_for_ccs_reclamation",
                    cmd = "filter",
                    args = '-tag "rq":"<' + min_read_quality + '" -length "<=' + max_reclamation_length + '"',
                    runtime_attr_override = disable_preemption_runtime_attrs
            }

            # Annotate our CCS uncorrected (reclaimable) reads
            call LONGBOW.Annotate as t_18_AnnotateReclaimableReads {
                input:
                    reads = t_17_ExtractCcsReclaimableReads.bam_out,
                    model = mas_seq_model,
                    runtime_attr_override = new_longbow_attrs,
            }

            call PB.PBIndex as t_19_PbIndexLongbowAnnotatedReclaimedReads {
                input:
                    bam = t_18_AnnotateReclaimableReads.annotated_bam
            }

            # 6: Longbow filter ccs reclaimable reads
            call LONGBOW.Filter as t_20_FilterReclaimableReads {
                input:
                    bam = t_18_AnnotateReclaimableReads.annotated_bam,
                    bam_pbi = t_19_PbIndexLongbowAnnotatedReclaimedReads.pbi,
                    prefix = sample_name + "_subshard",
                    model = mas_seq_model,
                    runtime_attr_override = new_longbow_attrs,
            }

            call PB.PBIndex as t_21_PbIndexLongbowAnnotatedReclaimedPassedReads {
                input:
                    bam = t_20_FilterReclaimableReads.passed_reads
            }

            # Shard these reads even wider so we can make sure we don't run out of memory:
            call PB.ShardLongReads as t_22_ShardReclaimedReads {
                input:
                    unaligned_bam = t_20_FilterReclaimableReads.passed_reads,
                    unaligned_pbi = t_21_PbIndexLongbowAnnotatedReclaimedPassedReads.pbi,
                    prefix = sample_name + "_ccs_reclaimed_longbow_annotated_subshard",
                    num_shards = secondary_scatter_width,
            }

            # Segment our arrays into individual array elements:
            scatter (corrected_shard in t_22_ShardReclaimedReads.unmapped_shards) {
                call LONGBOW.Segment as t_23_SegmentReclaimedAnnotatedReads {
                    input:
                        annotated_reads = corrected_shard,
                        model = mas_seq_model,
                        runtime_attr_override = new_longbow_attrs,
                        extra_args = segment_extra_args
                }
            }

            # Merge all outputs of Longbow Annotate / Segment:
            call Utils.MergeBams as t_24_MergeReclaimedArrayElements_1 {
                input:
                    bams = t_23_SegmentReclaimedAnnotatedReads.segmented_bam,
                    prefix = sample_name + "_ccs_reclaimed_ArrayElements_intermediate_1"
            }

            #######################################################################
            ## Align the data:
            #######################################################################

            # Align our array elements to the transcriptome:
            call PB.Align as t_12_AlignCCSArrayElementsToTranscriptome {
                input:
                    bam         = t_11_MergeCCSArrayElements_1.merged_bam,
                    ref_fasta   = transcriptome_ref_map['fasta'],
                    sample_name = sample_name,
                    map_preset  = "CCS"
            }

            # Align our array elements to the genome in splice-aware / ISOSEQ mode:
            call PB.Align as t_13_AlignCCSArrayElementsToGenome {
                input:
                    bam         = t_11_MergeCCSArrayElements_1.merged_bam,
                    ref_fasta   = genome_ref_map['fasta'],
                    sample_name = sample_name,
                    map_preset  = "ISOSEQ"
            }

            # To properly count our transcripts we must throw away the non-primary and unaligned reads:
            # Filter out all non-primary transcriptome alignments:
            call Utils.FilterReadsBySamFlags as t_14_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads {
                input:
                    bam = t_12_AlignCCSArrayElementsToTranscriptome.aligned_bam,
                    sam_flags = "2308",
                    prefix = sample_name + "_ccs_corrected_ArrayElements_Annotated_Transcriptome_Aligned_PrimaryOnly",
                    runtime_attr_override = filterReadsAttrs
            }
            # Filter out all non-primary genome alignments:
            call Utils.FilterReadsBySamFlags as t_15_RemoveUnmappedAndNonPrimaryCCSGenomeReads {
                input:
                    bam = t_13_AlignCCSArrayElementsToGenome.aligned_bam,
                    sam_flags = "2308",
                    prefix = sample_name + "_ccs_corrected_ArrayElements_Annotated_Genome_Aligned_PrimaryOnly",
                    runtime_attr_override = filterReadsAttrs
            }

            # Align our array elements to the transcriptome:
            # NOTE: We use the SUBREAD alignment preset here because these are uncorrected.
            call PB.Align as t_25_AlignReclaimedArrayElementsToTranscriptome {
                input:
                    bam         = t_24_MergeReclaimedArrayElements_1.merged_bam,
                    ref_fasta   = transcriptome_ref_map['fasta'],
                    sample_name = sample_name,
                    map_preset  = "SUBREAD"
            }

            # Align our array elements to the genome in splice-aware / ISOSEQ mode:
            # NOTE: We use the SUBREAD alignment preset here because these are uncorrected and we
            #       augment this with the exact parameters for ISOSEQ mode excapt we don't compress Homopolymer repeats:
            call PB.Align as t_26_AlignReclaimedArrayElementsToGenome {
                input:
                    bam         = t_24_MergeReclaimedArrayElements_1.merged_bam,
                    ref_fasta   = genome_ref_map['fasta'],
                    sample_name = sample_name,
                    map_preset  = "SUBREAD",
                    extra_options = "-k 15 -w 5 -o 2 -O 32 -e 1 -E 0 -A 1 -B 2 -z 200 -Z 100 -r 200000 -L 0.5 -g 2000 -C 5 -G 200000",
            }

            # To properly count our transcripts we must throw away the non-primary and unaligned reads:
            # Filter out all non-primary transcriptome alignments:
            call Utils.FilterReadsBySamFlags as t_27_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads {
                input:
                    bam = t_25_AlignReclaimedArrayElementsToTranscriptome.aligned_bam,
                    sam_flags = "2308",
                    prefix = sample_name + "_ccs_reclaimed_ArrayElements_Annotated_Transcriptome_Aligned_PrimaryOnly",
                    runtime_attr_override = filterReadsAttrs
            }
            # Filter out all non-primary genome alignments:
            call Utils.FilterReadsBySamFlags as t_28_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads {
                input:
                    bam = t_26_AlignReclaimedArrayElementsToGenome.aligned_bam,
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

        #################################################
        # CCS Passed:
        call Utils.MergeBams as t_29_MergeCCSReads { input: bams = t_05_FilterByMinQual.bam_out, prefix = sample_name + "_ccs_reads" }
        call Utils.MergeBams as t_30_MergeAnnotatedCCSReads { input: bams = t_06_LongbowAnnotateCCSReads.annotated_bam, prefix = sample_name + "_ccs_reads_annotated" }

        # Merge all CCS bams together for this Subread BAM:
        call Utils.MergeBams as t_31_MergePassedCCSReads { input: bams = t_07_FilterCCSReads.passed_reads, prefix = sample_name + "_ccs_longbow_passed_reads" }
        call Utils.MergeBams as t_32_MergeFailedCCSReads { input: bams = t_07_FilterCCSReads.failed_reads, prefix = sample_name + "_ccs_longbow_failed_reads" }
        call Utils.MergeBams as t_33_MergeCCSArrayElements { input: bams = t_11_MergeCCSArrayElements_1.merged_bam, prefix = sample_name + "_ccs_array_elements" }
        call Utils.MergeBams as t_34_MergeTranscriptomeAlignedCCSArrayElements { input: bams = t_12_AlignCCSArrayElementsToTranscriptome.aligned_bam, prefix = sample_name + "_ccs_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_35_MergeGenomeAlignedCCSArrayElements { input: bams = t_13_AlignCCSArrayElementsToGenome.aligned_bam, prefix = sample_name + "_ccs_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_36_MergePrimaryTranscriptomeAlignedCCSArrayElements { input: bams = t_14_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads.output_bam, prefix = sample_name + "_ccs_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_37_MergePrimaryGenomeAlignedCCSArrayElements { input: bams = t_15_RemoveUnmappedAndNonPrimaryCCSGenomeReads.output_bam, prefix = sample_name + "_ccs_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

        #################################################
        # CCS Reclaimed:
        call Utils.MergeBams as t_38_MergeCCSReclaimableReads { input: bams = t_17_ExtractCcsReclaimableReads.bam_out, prefix = sample_name + "_reclaimable_reads" }
        call Utils.MergeBams as t_39_MergeAnnotatedCCSReclaimableReads { input: bams = t_18_AnnotateReclaimableReads.annotated_bam, prefix = sample_name + "_reclaimable_reads_annotated" }

        # Merge all CCS bams together for this Subread BAM:
        call Utils.MergeBams as t_40_MergePassedCCSReclaimedReads { input: bams = t_20_FilterReclaimableReads.passed_reads, prefix = sample_name + "_reclaimable_longbow_passed_reads" }
        call Utils.MergeBams as t_41_MergeFailedCCSReclaimableReads { input: bams = t_20_FilterReclaimableReads.failed_reads, prefix = sample_name + "_reclaimed_longbow_failed_reads" }
        call Utils.MergeBams as t_42_MergeCCSReclaimedArrayElements { input: bams = t_24_MergeReclaimedArrayElements_1.merged_bam, prefix = sample_name + "_reclaimed_array_elements" }
        call Utils.MergeBams as t_43_MergeTranscriptomeAlignedCCSReclaimedArrayElements { input: bams = t_25_AlignReclaimedArrayElementsToTranscriptome.aligned_bam, prefix = sample_name + "_reclaimed_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_44_MergeGenomeAlignedCCSReclaimedArrayElements { input: bams = t_26_AlignReclaimedArrayElementsToGenome.aligned_bam, prefix = sample_name + "_reclaimed_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_45_MergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements { input: bams = t_27_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads.output_bam, prefix = sample_name + "_reclaimed_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_46_MergePrimaryGenomeAlignedCCSReclaimedArrayElements { input: bams = t_28_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads.output_bam, prefix = sample_name + "_reclaimed_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

        #################################################
        # Overall:
        call Utils.MergeBams as t_47_MergeAllAnnotatedReads { input: bams = flatten([t_06_LongbowAnnotateCCSReads.annotated_bam, t_18_AnnotateReclaimableReads.annotated_bam]), prefix = sample_name + "_all_annotated_reads" }
        call Utils.MergeBams as t_48_MergeAllLongbowPassedReads { input: bams = flatten([t_07_FilterCCSReads.passed_reads, t_20_FilterReclaimableReads.passed_reads]), prefix = sample_name + "_all_longbow_passed_reads" }
        call Utils.MergeBams as t_49_MergeAllLongbowFailedReads { input: bams = flatten([t_07_FilterCCSReads.failed_reads, t_20_FilterReclaimableReads.failed_reads]), prefix = sample_name + "_all_longbow_failed_reads" }

        call Utils.MergeBams as t_50_MergeAllArrayElements { input: bams = flatten([t_11_MergeCCSArrayElements_1.merged_bam, t_24_MergeReclaimedArrayElements_1.merged_bam]), prefix = sample_name + "_all_array_elements", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_51_MergeAllTranscriptomeAlignedArrayElements { input: bams = flatten([t_12_AlignCCSArrayElementsToTranscriptome.aligned_bam, t_25_AlignReclaimedArrayElementsToTranscriptome.aligned_bam]), prefix = sample_name + "_all_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_52_MergeAllGenomeAlignedArrayElements { input: bams = flatten([t_13_AlignCCSArrayElementsToGenome.aligned_bam, t_26_AlignReclaimedArrayElementsToGenome.aligned_bam]), prefix = sample_name + "_all_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_53_MergeAllPrimaryTranscriptomeAlignedArrayElements { input: bams = flatten([t_14_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads.output_bam, t_27_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads.output_bam]), prefix = sample_name + "_all_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_54_MergeAllPrimaryGenomeAlignedArrayElements { input: bams = flatten([t_15_RemoveUnmappedAndNonPrimaryCCSGenomeReads.output_bam, t_28_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads.output_bam]), prefix = sample_name + "_all_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    }

    #################################################
    #   _____ _             _
    #  |  ___(_)_ __   __ _| |
    #  | |_  | | '_ \ / _` | |
    #  |  _| | | | | | (_| | |__
    #  |_|   |_|_| |_|\__,_|____|
    #
    #   __  __
    #  |  \/  | ___ _ __ __ _  ___
    #  | |\/| |/ _ \ '__/ _` |/ _ \
    #  | |  | |  __/ | | (_| |  __/
    #  |_|  |_|\___|_|  \__, |\___|
    #                   |___/
    #
    # Now we have to merge all the files from each
    # input bam file together into one final output.
    #################################################

    #################################################
    # CCS Passed:
    call Utils.MergeBams as t_55_MasterMergeCCSReads { input: bams = t_29_MergeCCSReads.merged_bam, prefix = sample_name + "_ccs_reads" }
    call Utils.MergeBams as t_56_MasterMergeAnnotatedCCSReads { input: bams = t_30_MergeAnnotatedCCSReads.merged_bam, prefix = sample_name + "_ccs_reads_annotated" }

    # Merge all CCS bams together for this Subread BAM:
    call Utils.MergeBams as t_57_MasterMergePassedCCSReads { input: bams = t_31_MergePassedCCSReads.merged_bam, prefix = sample_name + "_ccs_longbow_passed_reads" }
    call Utils.MergeBams as t_58_MasterMergeFailedCCSReads { input: bams = t_32_MergeFailedCCSReads.merged_bam, prefix = sample_name + "_ccs_longbow_failed_reads" }
    call Utils.MergeBams as t_59_MasterMergeCCSArrayElements { input: bams = t_33_MergeCCSArrayElements.merged_bam, prefix = sample_name + "_ccs_array_elements" }
    call Utils.MergeBams as t_60_MasterMergeTranscriptomeAlignedCCSArrayElements { input: bams = t_34_MergeTranscriptomeAlignedCCSArrayElements.merged_bam, prefix = sample_name + "_ccs_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_61_MasterMergeGenomeAlignedCCSArrayElements { input: bams = t_35_MergeGenomeAlignedCCSArrayElements.merged_bam, prefix = sample_name + "_ccs_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_62_MasterMergePrimaryTranscriptomeAlignedCCSArrayElements { input: bams = t_36_MergePrimaryTranscriptomeAlignedCCSArrayElements.merged_bam, prefix = sample_name + "_ccs_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_63_MasterMergePrimaryGenomeAlignedCCSArrayElements { input: bams = t_37_MergePrimaryGenomeAlignedCCSArrayElements.merged_bam, prefix = sample_name + "_ccs_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    #################################################
    # CCS Reclaimed:
    call Utils.MergeBams as t_64_MasterMergeCCSReclaimableReads { input: bams = t_38_MergeCCSReclaimableReads.merged_bam, prefix = sample_name + "_reclaimable_reads" }
    call Utils.MergeBams as t_65_MasterMergeAnnotatedCCSReclaimableReads { input: bams = t_39_MergeAnnotatedCCSReclaimableReads.merged_bam, prefix = sample_name + "_reclaimable_reads_annotated" }

    # Merge all CCS bams together for this Subread BAM:
    call Utils.MergeBams as t_66_MasterMergePassedCCSReclaimedReads { input: bams = t_40_MergePassedCCSReclaimedReads.merged_bam, prefix = sample_name + "_reclaimable_longbow_passed_reads" }
    call Utils.MergeBams as t_67_MasterMergeFailedCCSReclaimableReads { input: bams = t_41_MergeFailedCCSReclaimableReads.merged_bam, prefix = sample_name + "_reclaimed_longbow_failed_reads" }
    call Utils.MergeBams as t_68_MasterMergeCCSReclaimedArrayElements { input: bams = t_42_MergeCCSReclaimedArrayElements.merged_bam, prefix = sample_name + "_reclaimed_array_elements" }
    call Utils.MergeBams as t_69_MasterMergeTranscriptomeAlignedCCSReclaimedArrayElements { input: bams = t_43_MergeTranscriptomeAlignedCCSReclaimedArrayElements.merged_bam, prefix = sample_name + "_reclaimed_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_70_MasterMergeGenomeAlignedCCSReclaimedArrayElements { input: bams = t_44_MergeGenomeAlignedCCSReclaimedArrayElements.merged_bam, prefix = sample_name + "_reclaimed_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_71_MasterMergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements { input: bams = t_45_MergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements.merged_bam, prefix = sample_name + "_reclaimed_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_72_MasterMergePrimaryGenomeAlignedCCSReclaimedArrayElements { input: bams = t_46_MergePrimaryGenomeAlignedCCSReclaimedArrayElements.merged_bam, prefix = sample_name + "_reclaimed_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    #################################################
    # Overall:
    call Utils.MergeBams as t_73_MasterMergeAllAnnotatedReads { input: bams = t_47_MergeAllAnnotatedReads.merged_bam, prefix = sample_name + "_all_annotated_reads" }
    call Utils.MergeBams as t_74_MasterMergeAllLongbowPassedReads { input: bams = t_48_MergeAllLongbowPassedReads.merged_bam, prefix = sample_name + "_all_longbow_passed_reads" }
    call Utils.MergeBams as t_75_MasterMergeAllLongbowFailedReads { input: bams = t_49_MergeAllLongbowFailedReads.merged_bam, prefix = sample_name + "_all_longbow_failed_reads" }

    call Utils.MergeBams as t_76_MasterMergeAllArrayElements { input: bams = t_50_MergeAllArrayElements.merged_bam, prefix = sample_name + "_all_array_elements", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_77_MasterMergeAllTranscriptomeAlignedArrayElements { input: bams = t_51_MergeAllTranscriptomeAlignedArrayElements.merged_bam, prefix = sample_name + "_all_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_78_MasterMergeAllGenomeAlignedArrayElements { input: bams = t_52_MergeAllGenomeAlignedArrayElements.merged_bam, prefix = sample_name + "_all_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_79_MasterMergeAllPrimaryTranscriptomeAlignedArrayElements { input: bams = t_53_MergeAllPrimaryTranscriptomeAlignedArrayElements.merged_bam, prefix = sample_name + "_all_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_80_MasterMergeAllPrimaryGenomeAlignedArrayElements { input: bams = t_54_MergeAllPrimaryGenomeAlignedArrayElements.merged_bam, prefix = sample_name + "_all_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    #################################################
    #         ___      ____
    #        / _ \    / ___|
    #       | | | |  | |
    #       | |_| |  | |___
    #        \__\_\   \____|
    #
    #################################################

    # Get stats on CCS reads:
    call LONGBOW.Stats as t_81_CCS_longbow_stats {
        input:
            reads = t_56_MasterMergeAnnotatedCCSReads.merged_bam,
            model = mas_seq_model,
            prefix = sample_name + "_CCS_Corrected",
            runtime_attr_override = new_longbow_attrs,
    }

    # Get stats on Reclaimable reads:
    call LONGBOW.Stats as t_82_Reclaimable_longbow_stats {
        input:
            reads = t_65_MasterMergeAnnotatedCCSReclaimableReads.merged_bam,
            model = mas_seq_model,
            prefix = sample_name + "_CCS_Reclaimable",
            runtime_attr_override = new_longbow_attrs,
    }

    # Get stats on Reclaimed reads:
    call LONGBOW.Stats as t_83_Reclaimed_longbow_stats {
        input:
            reads = t_66_MasterMergePassedCCSReclaimedReads.merged_bam,
            model = mas_seq_model,
            prefix = sample_name + "_CCS_Reclaimed",
            runtime_attr_override = new_longbow_attrs,
    }

    # Get stats on All reads (overall stats):
    call LONGBOW.Stats as t_84_Passed_longbow_stats {
        input:
            reads = t_74_MasterMergeAllLongbowPassedReads.merged_bam,
            model = mas_seq_model,
            prefix = sample_name + "_All_Longbow_Passed",
            runtime_attr_override = new_longbow_attrs,
    }

    # Get stats on All reads (overall stats):
    call LONGBOW.Stats as t_85_Overall_longbow_stats {
        input:
            reads = t_73_MasterMergeAllAnnotatedReads.merged_bam,
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

    String array_element_dir = outdir + "/array_elements"
    String arrays_dir = outdir + "/array_reads"

    ##############################################################################################################
    # Finalize the CCS Array Reads
    call FF.FinalizeToDir as t_86_FinalizeCCSReads {
        input:
            files = [
                t_55_MasterMergeCCSReads.merged_bam,
                t_55_MasterMergeCCSReads.merged_bai,
                t_56_MasterMergeAnnotatedCCSReads.merged_bam,
                t_56_MasterMergeAnnotatedCCSReads.merged_bai,
            ],
            outdir = arrays_dir,
            keyfile = t_85_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Reclaimed Array Reads
    call FF.FinalizeToDir as t_87_FinalizeCCSReclaimedReads {
        input:
            files = [
                t_64_MasterMergeCCSReclaimableReads.merged_bam,
                t_64_MasterMergeCCSReclaimableReads.merged_bai,
                t_65_MasterMergeAnnotatedCCSReclaimableReads.merged_bam,
                t_65_MasterMergeAnnotatedCCSReclaimableReads.merged_bai,
                t_66_MasterMergePassedCCSReclaimedReads.merged_bam,
                t_66_MasterMergePassedCCSReclaimedReads.merged_bai,
                t_67_MasterMergeFailedCCSReclaimableReads.merged_bam,
                t_67_MasterMergeFailedCCSReclaimableReads.merged_bai,
            ],
            outdir = arrays_dir,
            keyfile = t_85_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Overall Array Reads
    call FF.FinalizeToDir as t_88_FinalizeOverallCombinedReads {
        input:
            files = [
                t_73_MasterMergeAllAnnotatedReads.merged_bam,
                t_73_MasterMergeAllAnnotatedReads.merged_bai,
                t_74_MasterMergeAllLongbowPassedReads.merged_bam,
                t_74_MasterMergeAllLongbowPassedReads.merged_bai,
                t_75_MasterMergeAllLongbowFailedReads.merged_bam,
                t_75_MasterMergeAllLongbowFailedReads.merged_bai,
            ],
            outdir = arrays_dir,
            keyfile = t_85_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the CCS Array Element files
    call FF.FinalizeToDir as t_89_FinalizeCCSArrayElements {
        input:
            files = [
                t_59_MasterMergeCCSArrayElements.merged_bam,
                t_59_MasterMergeCCSArrayElements.merged_bai,
                t_60_MasterMergeTranscriptomeAlignedCCSArrayElements.merged_bam,
                t_60_MasterMergeTranscriptomeAlignedCCSArrayElements.merged_bai,
                t_61_MasterMergeGenomeAlignedCCSArrayElements.merged_bam,
                t_61_MasterMergeGenomeAlignedCCSArrayElements.merged_bai,
                t_62_MasterMergePrimaryTranscriptomeAlignedCCSArrayElements.merged_bam,
                t_62_MasterMergePrimaryTranscriptomeAlignedCCSArrayElements.merged_bai,
                t_63_MasterMergePrimaryGenomeAlignedCCSArrayElements.merged_bam,
                t_63_MasterMergePrimaryGenomeAlignedCCSArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = t_85_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Reclaimed Array Element files
    call FF.FinalizeToDir as t_90_FinalizeCCSRecaimedArrayElements {
        input:
            files = [
                t_68_MasterMergeCCSReclaimedArrayElements.merged_bam,
                t_68_MasterMergeCCSReclaimedArrayElements.merged_bai,
                t_69_MasterMergeTranscriptomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_69_MasterMergeTranscriptomeAlignedCCSReclaimedArrayElements.merged_bai,
                t_70_MasterMergeGenomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_70_MasterMergeGenomeAlignedCCSReclaimedArrayElements.merged_bai,
                t_71_MasterMergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_71_MasterMergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements.merged_bai,
                t_72_MasterMergePrimaryGenomeAlignedCCSReclaimedArrayElements.merged_bam,
                t_72_MasterMergePrimaryGenomeAlignedCCSReclaimedArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = t_85_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the Overall Array Element files
    call FF.FinalizeToDir as t_91_FinalizeOverallArrayElements {
        input:
            files = [
                t_76_MasterMergeAllArrayElements.merged_bam,
                t_76_MasterMergeAllArrayElements.merged_bai,
                t_77_MasterMergeAllTranscriptomeAlignedArrayElements.merged_bam,
                t_77_MasterMergeAllTranscriptomeAlignedArrayElements.merged_bai,
                t_78_MasterMergeAllGenomeAlignedArrayElements.merged_bam,
                t_78_MasterMergeAllGenomeAlignedArrayElements.merged_bai,
                t_79_MasterMergeAllPrimaryTranscriptomeAlignedArrayElements.merged_bam,
                t_79_MasterMergeAllPrimaryTranscriptomeAlignedArrayElements.merged_bai,
                t_80_MasterMergeAllPrimaryGenomeAlignedArrayElements.merged_bam,
                t_80_MasterMergeAllPrimaryGenomeAlignedArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = t_85_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Finalize the different stats files to separate directories
    call FF.FinalizeToDir as t_92_FinalizeCCSLongbowStats {
        input:
            files = [
                t_81_CCS_longbow_stats.summary_stats,
                t_81_CCS_longbow_stats.array_length_counts_plot_png,
                t_81_CCS_longbow_stats.array_length_counts_plot_svg,
                t_81_CCS_longbow_stats.ligation_heatmap_nn_png,
                t_81_CCS_longbow_stats.ligation_heatmap_nn_svg,
                t_81_CCS_longbow_stats.ligation_heatmap_png,
                t_81_CCS_longbow_stats.ligation_heatmap_svg,
                t_81_CCS_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_81_CCS_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_81_CCS_longbow_stats.ligation_heatmap_reduced_png,
                t_81_CCS_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = outdir + "/longbow_stats/CCS_Corrected/",
            keyfile = t_85_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_93_FinalizeReclaimableLongbowStats {
        input:
            files = [
                t_82_Reclaimable_longbow_stats.summary_stats,
                t_82_Reclaimable_longbow_stats.array_length_counts_plot_png,
                t_82_Reclaimable_longbow_stats.array_length_counts_plot_svg,
                t_82_Reclaimable_longbow_stats.ligation_heatmap_nn_png,
                t_82_Reclaimable_longbow_stats.ligation_heatmap_nn_svg,
                t_82_Reclaimable_longbow_stats.ligation_heatmap_png,
                t_82_Reclaimable_longbow_stats.ligation_heatmap_svg,
                t_82_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_82_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_82_Reclaimable_longbow_stats.ligation_heatmap_reduced_png,
                t_82_Reclaimable_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = outdir + "/longbow_stats/CCS_Reclaimable/",
            keyfile = t_85_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_94_FinalizeReclaimedLongbowStats {
        input:
            files = [
                t_83_Reclaimed_longbow_stats.summary_stats,
                t_83_Reclaimed_longbow_stats.array_length_counts_plot_png,
                t_83_Reclaimed_longbow_stats.array_length_counts_plot_svg,
                t_83_Reclaimed_longbow_stats.ligation_heatmap_nn_png,
                t_83_Reclaimed_longbow_stats.ligation_heatmap_nn_svg,
                t_83_Reclaimed_longbow_stats.ligation_heatmap_png,
                t_83_Reclaimed_longbow_stats.ligation_heatmap_svg,
                t_83_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_83_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_83_Reclaimed_longbow_stats.ligation_heatmap_reduced_png,
                t_83_Reclaimed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = outdir + "/longbow_stats/CCS_Reclaimed/",
            keyfile = t_85_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_95_FinalizeOverallLongbowStats {
        input:
            files = [
                t_85_Overall_longbow_stats.summary_stats,
                t_85_Overall_longbow_stats.array_length_counts_plot_png,
                t_85_Overall_longbow_stats.array_length_counts_plot_svg,
                t_85_Overall_longbow_stats.ligation_heatmap_nn_png,
                t_85_Overall_longbow_stats.ligation_heatmap_nn_svg,
                t_85_Overall_longbow_stats.ligation_heatmap_png,
                t_85_Overall_longbow_stats.ligation_heatmap_svg,
                t_85_Overall_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_85_Overall_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_85_Overall_longbow_stats.ligation_heatmap_reduced_png,
                t_85_Overall_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = outdir + "/longbow_stats/Overall/",
            keyfile = t_85_Overall_longbow_stats.summary_stats
    }
    call FF.FinalizeToDir as t_96_FinalizeAllPassedLongbowStats {
        input:
            files = [
                t_84_Passed_longbow_stats.summary_stats,
                t_84_Passed_longbow_stats.array_length_counts_plot_png,
                t_84_Passed_longbow_stats.array_length_counts_plot_svg,
                t_84_Passed_longbow_stats.ligation_heatmap_nn_png,
                t_84_Passed_longbow_stats.ligation_heatmap_nn_svg,
                t_84_Passed_longbow_stats.ligation_heatmap_png,
                t_84_Passed_longbow_stats.ligation_heatmap_svg,
                t_84_Passed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_84_Passed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_84_Passed_longbow_stats.ligation_heatmap_reduced_png,
                t_84_Passed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = outdir + "/longbow_stats/All_Longbow_Passed/",
            keyfile = t_85_Overall_longbow_stats.summary_stats
    }

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good
    call FF.WriteCompletionFile as t_97_WriteCompletionFile {
        input:
            outdir = outdir + "/",
            keyfile =  t_85_Overall_longbow_stats.summary_stats
    }

    ######################################################################
    #
    #          ___        _               _
    #         / _ \ _   _| |_ _ __  _   _| |_ ___
    #        | | | | | | | __| '_ \| | | | __/ __|
    #        | |_| | |_| | |_| |_) | |_| | |_\__ \
    #         \___/ \__,_|\__| .__/ \__,_|\__|___/
    #                        |_|
    #
    #######################################################################

    output {
#        File ccs_bam = FinalizeBam.gcs_path
#        File ccs_pbi = FinalizePbi.gcs_path
#
#        File aligned_bam = FinalizeAlignedBam.gcs_path
#        File aligned_bai = FinalizeAlignedBai.gcs_path
#        File aligned_pbi = FinalizeAlignedPbi.gcs_path
    }
}
