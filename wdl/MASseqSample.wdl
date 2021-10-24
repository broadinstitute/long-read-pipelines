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

        String sample_name

        String mas_seq_model = "mas15"

        Boolean is_SIRV_data = false
        Boolean generate_count_matrix = false

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

    # Setup our reference information:
    Map[String, String] genome_ref_map = read_map(genome_ref_map_file)
    Map[String, String] transcriptome_ref_map = read_map(transcriptome_ref_map_file)

    # Make sure we have a good output directory:
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/MASseqSample/~{sample_name}"

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

    # Do some error checking before we begin:
    Int num_ccs_bams = length(ccs_bams)
    Int num_ccs_pbis = length(ccs_pbis)
    if (num_ccs_bams != num_ccs_pbis) {
        call Utils.FailWithWarning as t_04_WARN1 { input: warning = "Error: Num CCS bams is not the same as the number of CCS pbindex files  (~{num_ccs_bams} != ~{num_ccs_pbis}).  Aborting." }
    }

    # Loop through all of our base input files and process them individually.
    # We'll merge everything together at the end.
    scatter ( i in range(length(ccs_bams))) {
        File ccs_bam = ccs_bams[i]
        File ccs_pbi = ccs_pbis[i]

        call PB.ShardLongReads as t_08_ShardLongReads {
            input:
                unaligned_bam = ccs_bam,
                unaligned_pbi = ccs_pbi,
                prefix = sample_name + "_shard",
                num_shards = 300,
                runtime_attr_override = fast_sharding_runtime_attrs,
        }

        scatter (sharded_reads in t_08_ShardLongReads.unmapped_shards) {

            #######################################################################
            ## Setup
            #######################################################################

            String fbmrq_prefix = basename(ccs_bam, ".bam")

            # Filter out the kinetics tags from PB files:
            call PB.RemoveKineticsTags as t_07_RemoveKineticsTags {
                input:
                    bam = ccs_bam,
                    prefix = sample_name + "_kinetics_removed"
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
                    model = mas_seq_model,
                    runtime_attr_override = new_longbow_attrs,
            }

            # 6: Longbow filter ccs reclaimable reads
            call LONGBOW.Filter as t_10_FilterCCSReads {
                input:
                    bam = t_09_LongbowAnnotateCCSReads.annotated_bam,
                    prefix = sample_name + "_ccs_corrected_subshard",
                    model = mas_seq_model,
                    runtime_attr_override = new_longbow_attrs,
            }

            call PB.PBIndex as t_11_PbIndexLongbowAnnotatedCCSPassedReads {
                input:
                    bam = t_10_FilterCCSReads.passed_reads
            }

            # Shard these reads even wider so we can make sure we don't run out of memory:
            call PB.ShardLongReads as t_12_ShardCorrectedReads {
                input:
                    unaligned_bam = t_10_FilterCCSReads.passed_reads,
                    unaligned_pbi = t_11_PbIndexLongbowAnnotatedCCSPassedReads.pbi,
                    prefix = sample_name + "_ccs_corrected_longbow_annotated_subshard",
                    num_shards = secondary_scatter_width,
            }

            # Segment our arrays into individual array elements:
            scatter (corrected_shard in t_12_ShardCorrectedReads.unmapped_shards) {
                call LONGBOW.Segment as t_13_SegmentCCSAnnotatedReads {
                    input:
                        annotated_reads = corrected_shard,
                        model = mas_seq_model,
                        runtime_attr_override = new_longbow_attrs,
                }
            }

            # Merge all outputs of Longbow Annotate / Segment:
            call Utils.MergeBams as t_14_MergeCCSArrayElements_1 {
                input:
                    bams = t_13_SegmentCCSAnnotatedReads.segmented_bam,
                    prefix = sample_name + "_ccs_corrected_ArrayElements_intermediate_1"
            }

            # Align our array elements to the transcriptome:
            call PB.Align as t_15_AlignCCSArrayElementsToTranscriptome {
                input:
                    bam         = t_14_MergeCCSArrayElements_1.merged_bam,
                    ref_fasta   = transcriptome_ref_map['fasta'],
                    sample_name = sample_name,
                    map_preset  = "CCS"
            }

            # Align our array elements to the genome in splice-aware / ISOSEQ mode:
            call PB.Align as t_16_AlignCCSArrayElementsToGenome {
                input:
                    bam         = t_14_MergeCCSArrayElements_1.merged_bam,
                    ref_fasta   = genome_ref_map['fasta'],
                    sample_name = sample_name,
                    map_preset  = "ISOSEQ"
            }

            # To properly count our transcripts we must throw away the non-primary and unaligned reads:
            # Filter out all non-primary transcriptome alignments:
            call Utils.FilterReadsBySamFlags as t_19_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads {
                input:
                    bam = t_15_AlignCCSArrayElementsToTranscriptome.aligned_bam,
                    sam_flags = "2308",
                    prefix = sample_name + "_ccs_corrected_ArrayElements_Annotated_Transcriptome_Aligned_PrimaryOnly",
                    runtime_attr_override = filterReadsAttrs
            }
            # Filter out all non-primary genome alignments:
            call Utils.FilterReadsBySamFlags as t_20_RemoveUnmappedAndNonPrimaryCCSGenomeReads {
                input:
                    bam = t_16_AlignCCSArrayElementsToGenome.aligned_bam,
                    sam_flags = "2308",
                    prefix = sample_name + "_ccs_corrected_ArrayElements_Annotated_Genome_Aligned_PrimaryOnly",
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
                    model = mas_seq_model,
                    runtime_attr_override = new_longbow_attrs,
            }

            call PB.PBIndex as t_24_PbIndexLongbowAnnotatedReclaimedReads {
                input:
                    bam = t_23_AnnotateReclaimableReads.annotated_bam
            }

            # 6: Longbow filter ccs reclaimable reads
            call LONGBOW.Filter as t_25_FilterReclaimableReads {
                input:
                    bam = t_23_AnnotateReclaimableReads.annotated_bam,
                    bam_pbi = t_24_PbIndexLongbowAnnotatedReclaimedReads.pbi,
                    prefix = sample_name + "_subshard",
                    model = mas_seq_model,
                    runtime_attr_override = new_longbow_attrs,
            }

            call PB.PBIndex as t_26_PbIndexLongbowAnnotatedReclaimedPassedReads {
                input:
                    bam = t_25_FilterReclaimableReads.passed_reads
            }

            # Shard these reads even wider so we can make sure we don't run out of memory:
            call PB.ShardLongReads as t_27_ShardReclaimedReads {
                input:
                    unaligned_bam = t_25_FilterReclaimableReads.passed_reads,
                    unaligned_pbi = t_26_PbIndexLongbowAnnotatedReclaimedPassedReads.pbi,
                    prefix = sample_name + "_ccs_reclaimed_longbow_annotated_subshard",
                    num_shards = secondary_scatter_width,
            }

            # Segment our arrays into individual array elements:
            scatter (corrected_shard in t_27_ShardReclaimedReads.unmapped_shards) {
                call LONGBOW.Segment as t_28_SegmentReclaimedAnnotatedReads {
                    input:
                        annotated_reads = corrected_shard,
                        model = mas_seq_model,
                        runtime_attr_override = new_longbow_attrs,
                }
            }

            # Merge all outputs of Longbow Annotate / Segment:
            call Utils.MergeBams as t_29_MergeReclaimedArrayElements_1 {
                input:
                    bams = t_28_SegmentReclaimedAnnotatedReads.segmented_bam,
                    prefix = sample_name + "_ccs_reclaimed_ArrayElements_intermediate_1"
            }

            # Align our array elements to the transcriptome:
            # NOTE: We use the SUBREAD alignment preset here because these are uncorrected.
            call PB.Align as t_30_AlignReclaimedArrayElementsToTranscriptome {
                input:
                    bam         = t_29_MergeReclaimedArrayElements_1.merged_bam,
                    ref_fasta   = transcriptome_ref_map['fasta'],
                    sample_name = sample_name,
                    map_preset  = "SUBREAD"
            }

            # Align our array elements to the genome in splice-aware / ISOSEQ mode:
            # NOTE: We use the SUBREAD alignment preset here because these are uncorrected and we
            #       augment this with the exact parameters for ISOSEQ mode excapt we don't compress Homopolymer repeats:
            call PB.Align as t_31_AlignReclaimedArrayElementsToGenome {
                input:
                    bam         = t_29_MergeReclaimedArrayElements_1.merged_bam,
                    ref_fasta   = genome_ref_map['fasta'],
                    sample_name = sample_name,
                    map_preset  = "SUBREAD",
                    extra_options = "-k 15 -w 5 -o 2 -O 32 -e 1 -E 0 -A 1 -B 2 -z 200 -Z 100 -r 200000 -L 0.5 -g 2000 -C 5 -G 200000",
            }

            # To properly count our transcripts we must throw away the non-primary and unaligned reads:
            # Filter out all non-primary transcriptome alignments:
            call Utils.FilterReadsBySamFlags as t_34_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads {
                input:
                    bam = t_30_AlignReclaimedArrayElementsToTranscriptome.aligned_bam,
                    sam_flags = "2308",
                    prefix = sample_name + "_ccs_reclaimed_ArrayElements_Annotated_Transcriptome_Aligned_PrimaryOnly",
                    runtime_attr_override = filterReadsAttrs
            }
            # Filter out all non-primary genome alignments:
            call Utils.FilterReadsBySamFlags as t_35_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads {
                input:
                    bam = t_31_AlignReclaimedArrayElementsToGenome.aligned_bam,
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
        call Utils.MergeBams as t_36_MergeCCSReads { input: bams = t_08_FilterByMinQual.bam_out, prefix = sample_name + "_ccs_reads" }
        call Utils.MergeBams as t_37_MergeAnnotatedCCSReads { input: bams = t_09_LongbowAnnotateCCSReads.annotated_bam, prefix = sample_name + "_ccs_reads_annotated" }

        # Merge all CCS bams together for this Subread BAM:
        call Utils.MergeBams as t_38_MergePassedCCSReads { input: bams = t_10_FilterCCSReads.passed_reads, prefix = sample_name + "_ccs_longbow_passed_reads" }
        call Utils.MergeBams as t_39_MergeFailedCCSReads { input: bams = t_10_FilterCCSReads.failed_reads, prefix = sample_name + "_ccs_longbow_failed_reads" }
        call Utils.MergeBams as t_40_MergeCCSArrayElements { input: bams = t_14_MergeCCSArrayElements_1.merged_bam, prefix = sample_name + "_ccs_array_elements" }
        call Utils.MergeBams as t_41_MergeTranscriptomeAlignedCCSArrayElements { input: bams = t_15_AlignCCSArrayElementsToTranscriptome.aligned_bam, prefix = sample_name + "_ccs_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_42_MergeGenomeAlignedCCSArrayElements { input: bams = t_16_AlignCCSArrayElementsToGenome.aligned_bam, prefix = sample_name + "_ccs_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_43_MergePrimaryTranscriptomeAlignedCCSArrayElements { input: bams = t_19_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads.output_bam, prefix = sample_name + "_ccs_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_44_MergePrimaryGenomeAlignedCCSArrayElements { input: bams = t_20_RemoveUnmappedAndNonPrimaryCCSGenomeReads.output_bam, prefix = sample_name + "_ccs_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

        #################################################
        # CCS Reclaimed:
        call Utils.MergeBams as t_45_MergeCCSReclaimableReads { input: bams = t_22_ExtractCcsReclaimableReads.bam_out, prefix = sample_name + "_reclaimable_reads" }
        call Utils.MergeBams as t_46_MergeAnnotatedCCSReclaimableReads { input: bams = t_23_AnnotateReclaimableReads.annotated_bam, prefix = sample_name + "_reclaimable_reads_annotated" }

        # Merge all CCS bams together for this Subread BAM:
        call Utils.MergeBams as t_47_MergePassedCCSReclaimedReads { input: bams = t_25_FilterReclaimableReads.passed_reads, prefix = sample_name + "_reclaimable_longbow_passed_reads" }
        call Utils.MergeBams as t_48_MergeFailedCCSReclaimableReads { input: bams = t_25_FilterReclaimableReads.failed_reads, prefix = sample_name + "_reclaimed_longbow_failed_reads" }
        call Utils.MergeBams as t_49_MergeCCSReclaimedArrayElements { input: bams = t_29_MergeReclaimedArrayElements_1.merged_bam, prefix = sample_name + "_reclaimed_array_elements" }
        call Utils.MergeBams as t_50_MergeTranscriptomeAlignedCCSReclaimedArrayElements { input: bams = t_30_AlignReclaimedArrayElementsToTranscriptome.aligned_bam, prefix = sample_name + "_reclaimed_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_51_MergeGenomeAlignedCCSReclaimedArrayElements { input: bams = t_31_AlignReclaimedArrayElementsToGenome.aligned_bam, prefix = sample_name + "_reclaimed_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_52_MergePrimaryTranscriptomeAlignedCCSReclaimedArrayElements { input: bams = t_34_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads.output_bam, prefix = sample_name + "_reclaimed_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_53_MergePrimaryGenomeAlignedCCSReclaimedArrayElements { input: bams = t_35_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads.output_bam, prefix = sample_name + "_reclaimed_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

        #################################################
        # Overall:
        call Utils.MergeBams as t_54_MergeAllAnnotatedReads { input: bams = flatten([t_09_LongbowAnnotateCCSReads.annotated_bam, t_23_AnnotateReclaimableReads.annotated_bam]), prefix = sample_name + "_all_annotated_reads" }
        call Utils.MergeBams as t_55_MergeAllLongbowPassedReads { input: bams = flatten([t_10_FilterCCSReads.passed_reads, t_25_FilterReclaimableReads.passed_reads]), prefix = sample_name + "_all_longbow_passed_reads" }
        call Utils.MergeBams as t_56_MergeAllLongbowFailedReads { input: bams = flatten([t_10_FilterCCSReads.failed_reads, t_25_FilterReclaimableReads.failed_reads]), prefix = sample_name + "_all_longbow_failed_reads" }

        call Utils.MergeBams as t_57_MergeAllArrayElements { input: bams = flatten([t_14_MergeCCSArrayElements_1.merged_bam, t_29_MergeReclaimedArrayElements_1.merged_bam]), prefix = sample_name + "_all_array_elements", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_58_MergeAllTranscriptomeAlignedArrayElements { input: bams = flatten([t_15_AlignCCSArrayElementsToTranscriptome.aligned_bam, t_30_AlignReclaimedArrayElementsToTranscriptome.aligned_bam]), prefix = sample_name + "_all_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_59_MergeAllGenomeAlignedArrayElements { input: bams = flatten([t_16_AlignCCSArrayElementsToGenome.aligned_bam, t_31_AlignReclaimedArrayElementsToGenome.aligned_bam]), prefix = sample_name + "_all_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_60_MergeAllPrimaryTranscriptomeAlignedArrayElements { input: bams = flatten([t_19_RemoveUnmappedAndNonPrimaryCCSTranscriptomeReads.output_bam, t_34_RemoveUnmappedAndNonPrimaryTranscriptomeReclaimedReads.output_bam]), prefix = sample_name + "_all_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
        call Utils.MergeBams as t_61_MergeAllPrimaryGenomeAlignedArrayElements { input: bams = flatten([t_20_RemoveUnmappedAndNonPrimaryCCSGenomeReads.output_bam, t_35_RemoveUnmappedAndNonPrimaryGenomeReclaimedReads.output_bam]), prefix = sample_name + "_all_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    }
}
