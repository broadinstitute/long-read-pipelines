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

workflow MasSeqArraySlideSeq {

    meta {
        description : "This workflow will process and demultiplex data from the MAS-seq Slide Seq array."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String gcs_input_dir
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/MasSeqArraySlideSeq"

        File head_adapter_fasta = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/10x_adapter.fasta"
        File tail_adapter_fasta = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/tso_adapter.fasta"
        File ten_x_cell_barcode_whitelist = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/10x_Barcodes_3M-february-2018.txt"

        # NOTE: Reference for un-split CCS reads:
        File ref_fasta =  "gs://broad-dsde-methods-long-reads/resources/references/grch38/Homo_sapiens_assembly38.fasta"
        File ref_fasta_index = "gs://broad-dsde-methods-long-reads/resources/references/grch38/Homo_sapiens_assembly38.fasta.fai"
        File ref_fasta_dict = "gs://broad-dsde-methods-long-reads/resources/references/grch38/Homo_sapiens_assembly38.dict"

        # NOTE: Reference for array elements:
        File transcriptome_ref_fasta =  "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.pc_transcripts.fa"
        File transcriptome_ref_fasta_index = "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.pc_transcripts.fa.fai"
        File transcriptome_ref_fasta_dict = "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.pc_transcripts.dict"

        File genome_annotation_gtf = "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.primary_assembly.annotation.gtf"

        File jupyter_template_static = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/MAS-seq_QC_report_template-static.ipynb"
        File workflow_dot_file = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/PB10xMasSeqArraySingleFlowcellv2.dot"

        # Default here is 0 because ccs uncorrected reads all seem to have RQ = -1.
        # All pathologically long reads also have RQ = -1.
        # This way we preserve the vast majority of the data, even if it has low quality.
        # We can filter it out at later steps.
        Float min_read_quality = 0.0
        Int max_reclamation_length = 60000

        String? sample_name
    }

    parameter_meta {
        gcs_input_dir : "Input folder on GCS in which to search for BAM files to process."
        gcs_out_root_dir : "Root output GCS folder in which to place results of this workflow."

        head_adapter_fasta : "FASTA file containing the sequence that each transcript should start with.  Typically this will be the 10x adapter sequence from the 10x library prep."
        tail_adapter_fasta : "FASTA file containing the sequence that each transcript should end with.  Typically this will be the Template Switch Oligo (TSO) sequence from the 10x library prep."
        ten_x_cell_barcode_whitelist : "Text file containing a whitelist of cell barcodes for the 10x library prep."

        ref_fasta : "FASTA file containing the reference sequence to which the input data should be aligned before splitting into array elements."
        ref_fasta_index : "FASTA index file for the given ref_fasta file."
        ref_fasta_dict : "Sequence dictionary file for the given ref_fasta file."

        transcriptome_ref_fasta : "FASTA file containing the reference sequence to which the array elements should be aligned."
        transcriptome_ref_fasta_index : "FASTA index file for the given transcriptome_ref_fasta file."
        transcriptome_ref_fasta_dict : "Sequence dictionary file for the given transcriptome_ref_fasta file."

        genome_annotation_gtf : "Gencode GTF file containing genome annotations for the organism under study (usually humans).  This must match the given reference version (usually hg38)."

        jupyter_template_static : "Jupyter notebook / ipynb file containing a template for the QC report which will contain static plots.  This should contain the same information as the jupyter_template_interactive file, but with static images."
        workflow_dot_file : "DOT file containing the representation of this WDL to be included in the QC reports.  This can be generated with womtool."

        min_read_quality : "[optional] Minimum read quality for reads to have to be included in our data (Default: 0.0)."
        max_reclamation_length : "[optional] Maximum length (in bases) that a read can be to attempt to reclaim from CCS rejection (Default: 60000)."

        sample_name : "[optional] The name of the sample to associate with the data in this workflow."
    }

    # Set our model here:
    String mas_seq_model = "slide-seq"

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as WdlExecutionStartTimestamp { input: }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = gcs_input_dir }

    # Check here if we found ccs bams or subread bams:
    Boolean use_subreads = FindBams.has_subreads

    # Make sure we have **EXACTLY** one bam file to run on:
    if (length(FindBams.ccs_bams) != 1) {
        call Utils.FailWithWarning { input: warning = "Error: Multiple BAM files found.  Cannot continue!" }
     }

    # Alias our bam file so we can work with it easier:
    File reads_bam = FindBams.ccs_bams[0]

    call PB.GetRunInfo { input: subread_bam = reads_bam }

    String SM  = select_first([sample_name, GetRunInfo.run_info["SM"]])
    String PL  = "PACBIO"
    String PU  = GetRunInfo.run_info["PU"]
    String DT  = GetRunInfo.run_info["DT"]
    String ID  = PU
    String DS  = GetRunInfo.run_info["DS"]
    String DIR = SM + "." + ID

    File read_pbi = sub(reads_bam, ".bam$", ".bam.pbi")
    call PB.ShardLongReads {
        input:
            unaligned_bam = reads_bam,
            unaligned_pbi = read_pbi,
            prefix = SM + "_shard",
            num_shards = 50,
    }

    scatter (sharded_reads in ShardLongReads.unmapped_shards) {

        ## No more preemption on this sharding - takes too long otherwise.
        RuntimeAttr disable_preemption_runtime_attrs = object {
            preemptible_tries: 0
        }

        String fbmrq_prefix = basename(sharded_reads, ".bam")

        # Filter out the kinetics tags from PB files:
        call PB.RemoveKineticsTags {
            input:
                bam = sharded_reads,
                prefix = SM + "_kinetics_removed"
        }

        # Handle setting up the things that we need for further processing of CCS-only reads:
        call PB.FindCCSReport {
            input:
                gcs_input_dir = gcs_input_dir
        }

        # 1 - filter the reads by the minimum read quality:
        call Utils.Bamtools as FilterByMinQual {
            input:
                bamfile = RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_good_reads",
                cmd = "filter",
                args = '-tag "rq":">=' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # Get reclaimable reads:
        call Utils.Bamtools as GetCcsReclaimableReads {
            input:
                bamfile = RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_reads_for_ccs_reclamation",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '" -length "<=' + max_reclamation_length + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # Annotate our reads:
        call LONGBOW.Annotate as LongbowAnnotateReads {
            input:
                reads = FilterByMinQual.bam_out,
                model = mas_seq_model
        }

        # Annotate our reclaimable reads:
        call LONGBOW.Annotate as AnnotateCcsReclaimableReads {
            input:
                reads = GetCcsReclaimableReads.bam_out,
                model = mas_seq_model
        }

        # Longbow filter our reads:
        call LONGBOW.Filter as FilterLongbowAnnotatedReads {
            input:
                bam = LongbowAnnotateReads.annotated_bam,
                prefix = SM + "_longbow_filtered",
                model = mas_seq_model
        }
        call PB.PBIndex as PbIndexLongbowFilteredReads {
            input:
                bam = FilterLongbowAnnotatedReads.passed_reads
        }

        # Longbow filter ccs reclaimable reads to get CCS Reclaimed reads:
        call LONGBOW.Filter as FilterCcsReclaimableLongbowAnnotatedReads {
            input:
                bam = AnnotateCcsReclaimableReads,
                prefix = SM + "_ccs_reclaimable_longbow_filtered",
                model = mas_seq_model
        }
        File ccs_reclaimed_reads = FilterCcsReclaimableLongbowAnnotatedReads.passed_reads
        call PB.PBIndex as PbIndexLongbowFilteredCcsReclaimedReads {
            input:
                bam = ccs_reclaimed_reads
        }

        # Shard the ccs reclaimed reads even wider so we can make sure we don't run out of memory:
        call PB.ShardLongReads as ShardReclaimedReads {
            input:
                unaligned_bam = ccs_reclaimed_reads,
                unaligned_pbi = PbIndexLongbowFilteredCcsReclaimedReads.pbindex,
                prefix = SM + "_ccs_reclaimed_longbow_annotated_filtered_subshard",
                num_shards = 10,
        }

        # Shard the ccs corrected reads even wider so we can make sure we don't run out of memory:
        call PB.ShardLongReads as ShardCorrectedReads {
            input:
                unaligned_bam = FilterLongbowAnnotatedReads.passed_reads,
                unaligned_pbi = PbIndexLongbowFilteredReads.pbindex,
                prefix = SM + "_longbow_annotated_filtered_subshard",
                num_shards = 10,
        }

        # Segment our CCS Corrected arrays into individual array elements and align them:
        scatter (corrected_shard in ShardCorrectedReads.unmapped_shards) {
            call LONGBOW.Segment as SegmentAnnotatedReads {
                input:
                    annotated_reads = corrected_shard,
                    model = mas_seq_model
            }

            call AR.Minimap2 as AlignArrayElementsToGenome {
                input:
                    reads      = [ SegmentAnnotatedReads.segmented_bam ],
                    ref_fasta  = ref_fasta,
                    map_preset = "splice:hq"
            }

            # Align our array elements to the transcriptome:
            call AR.Minimap2 as AlignArrayElements {
                input:
                    reads      = [ SegmentAnnotatedReads.segmented_bam ],
                    ref_fasta  = transcriptome_ref_fasta,
                    map_preset = "splice:hq"
            }

            # We need to restore the annotations we created with the 10x tool to the aligned reads.
            call TENX.RestoreAnnotationstoAlignedBam as RestoreAnnotationsToTranscriptomeAlignedBam {
                input:
                    annotated_bam_file = SegmentAnnotatedReads.segmented_bam,
                    aligned_bam_file = AlignArrayElements.aligned_bam,
                    tags_to_ignore = [],
                    mem_gb = 32,  # TODO: Debug for memory redution (was 64)
            }

            # We need to restore the annotations we created with the 10x tool to the aligned reads.
            call TENX.RestoreAnnotationstoAlignedBam as RestoreAnnotationsToGenomeAlignedBam {
                input:
                    annotated_bam_file = SegmentAnnotatedReads.segmented_bam,
                    aligned_bam_file = AlignArrayElementsToGenome.aligned_bam,
                    tags_to_ignore = [],
                    mem_gb = 32,  # TODO: Debug for memory redution (was 64)
            }
        }

        # Merge all outputs of Longbow Annotate / Segment:
        call Utils.MergeBams as MergeRawArrayElementSubshards {
            input:
                bams = SegmentAnnotatedReads.segmented_bam,
                prefix = SM + "_ArrayElements_Filtered_intermediate_1"
        }
        call Utils.MergeBams as MergeTxAlignedArrayElementSubshards {
            input:
                bams = RestoreAnnotationsToTranscriptomeAlignedBam.output_bam,
                prefix = SM + "_ArrayElements_Annotated_Filtered_Transcriptome_Aligned"
        }
        call Utils.MergeBams as MergeGenomeAlignedArrayElementSubshards {
            input:
                bams = RestoreAnnotationsToGenomeAlignedBam.output_bam,
                prefix = SM + "_ArrayElements_Annotated_Filtered_Genome_Aligned"
        }

        # To properly count our transcripts we must throw away the non-primary and unaligned reads:
        RuntimeAttr filterReadsAttrs = object {
            cpu_cores: 4,
            preemptible_tries: 0
        }
        call Utils.FilterReadsBySamFlags as RemoveUnmappedAndNonPrimaryTranscriptomeReads {
            input:
                bam = MergeTxAlignedArrayElementSubshards.merged_bam,
                sam_flags = "2308",
                prefix = SM + "_ArrayElements_Annotated_Filtered_Transcriptome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }
        call Utils.FilterReadsBySamFlags as RemoveUnmappedAndNonPrimaryGenomeReads {
            input:
                bam = MergeGenomeAlignedArrayElementSubshards.merged_bam,
                sam_flags = "2308",
                prefix = SM + "_ArrayElements_Annotated_Filtered_Genome_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }

        # Segment our CCS Corrected arrays into individual array elements and align them:
        scatter (corrected_shard in ShardCorrectedReads.unmapped_shards) {
            call LONGBOW.Segment as SegmentAnnotatedReads {
                input:
                    annotated_reads = corrected_shard,
                    model = mas_seq_model
            }

            call AR.Minimap2 as AlignArrayElementsToGenome {
                input:
                    reads      = [ SegmentAnnotatedReads.segmented_bam ],
                    ref_fasta  = ref_fasta,
                    map_preset = "splice:hq"
            }

            # Align our array elements to the transcriptome:
            call AR.Minimap2 as AlignArrayElements {
                input:
                    reads      = [ SegmentAnnotatedReads.segmented_bam ],
                    ref_fasta  = transcriptome_ref_fasta,
                    map_preset = "splice:hq"
            }

            # We need to restore the annotations we created with the 10x tool to the aligned reads.
            call TENX.RestoreAnnotationstoAlignedBam as RestoreAnnotationsToTranscriptomeAlignedBam {
                input:
                    annotated_bam_file = SegmentAnnotatedReads.segmented_bam,
                    aligned_bam_file = AlignArrayElements.aligned_bam,
                    tags_to_ignore = [],
                    mem_gb = 32,  # TODO: Debug for memory redution (was 64)
            }

            # We need to restore the annotations we created with the 10x tool to the aligned reads.
            call TENX.RestoreAnnotationstoAlignedBam as RestoreAnnotationsToGenomeAlignedBam {
                input:
                    annotated_bam_file = SegmentAnnotatedReads.segmented_bam,
                    aligned_bam_file = AlignArrayElementsToGenome.aligned_bam,
                    tags_to_ignore = [],
                    mem_gb = 32,  # TODO: Debug for memory redution (was 64)
            }
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
    call Utils.MergeBams as MergeCCSRqFilteredReads { input: bams = FilterByMinQual.bam_out, prefix = SM + "_ccs_reads" }
    call Utils.MergeBams as MergeAnnotatedCCSReads { input: bams = LongbowAnnotateReads.annotated_bam, prefix = SM + "_ccs_reads_annotated" }
    call Utils.MergeBams as MergeLongbowFilteredPassedReads { input: bams = FilterLongbowAnnotatedReads.passed_reads, prefix = SM + "_ccs_reads_annotated_filtered_passed" }
    call Utils.MergeBams as MergeLongbowFilteredFailedReads { input: bams = FilterLongbowAnnotatedReads.failed_reads, prefix = SM + "_ccs_reads_annotated_filtered_failed" }

    # Merge all CCS bams together for this Subread BAM:
    RuntimeAttr merge_extra_cpu_attrs = object {
        cpu_cores: 4
    }
    call Utils.MergeBams as MergeRawArrayElements { input: bams = MergeRawArrayElementSubshards.merged_bam, prefix = SM + "_raw_array_elements" }
    call Utils.MergeBams as MergeTranscriptomeAlignedArrayElements { input: bams = MergeTxAlignedArrayElementSubshards.merged_bam, prefix = SM + "_array_elements_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as MergeGenomeAlignedArrayElements { input: bams = MergeGenomeAlignedArrayElementSubshards.merged_bam, prefix = SM + "_array_elements_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as MergePrimaryTranscriptomeAlignedArrayElements { input: bams = RemoveUnmappedAndNonPrimaryTranscriptomeReads.output_bam, prefix = SM + "_array_elements_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as MergePrimaryGenomeAlignedArrayElements { input: bams = RemoveUnmappedAndNonPrimaryGenomeReads.output_bam, prefix = SM + "_array_elements_genome_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

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
    String base_out_dir = outdir + "/" + DIR + "/" + WdlExecutionStartTimestamp.timestamp_string

    File final_ccs_report = FindCCSReport.ccs_report[0]

    String array_element_dir = base_out_dir + "/array_elements"
    String arrays_dir = base_out_dir + "/array_reads"

    ##############################################################################################################
    # Finalize the final annotated, aligned array elements:
    call FF.FinalizeToDir as FinalizeArrayReads {
        input:
            files = [
                MergeCCSRqFilteredReads.merged_bam,
                MergeCCSRqFilteredReads.merged_bai,
                MergeAnnotatedCCSReads.merged_bam,
                MergeAnnotatedCCSReads.merged_bai,
                MergeLongbowFilteredPassedReads.merged_bam,
                MergeLongbowFilteredPassedReads.merged_bai,
                MergeLongbowFilteredFailedReads.merged_bam,
                MergeLongbowFilteredFailedReads.merged_bai,
            ],
            outdir = arrays_dir,
            keyfile = MergePrimaryGenomeAlignedArrayElements.merged_bai
    }

    ##############################################################################################################
    # Finalize the intermediate reads files (from raw CCS corrected reads through split array elements)
    call FF.FinalizeToDir as FinalizeArrayElements {
        input:
            files = [
                MergeRawArrayElements.merged_bam,
                MergeRawArrayElements.merged_bai,
                MergeTranscriptomeAlignedArrayElements.merged_bam,
                MergeTranscriptomeAlignedArrayElements.merged_bai,
                MergeGenomeAlignedArrayElements.merged_bam,
                MergeGenomeAlignedArrayElements.merged_bai,
                MergePrimaryTranscriptomeAlignedArrayElements.merged_bam,
                MergePrimaryTranscriptomeAlignedArrayElements.merged_bai,
                MergePrimaryGenomeAlignedArrayElements.merged_bam,
                MergePrimaryGenomeAlignedArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = MergePrimaryGenomeAlignedArrayElements.merged_bai
    }

    call FF.FinalizeToDir as FinalizeCCSReport {
        input:
            files = [
                final_ccs_report
            ],
            outdir = base_out_dir + "/",
            keyfile = MergePrimaryGenomeAlignedArrayElements.merged_bai
    }

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good:
    call FF.WriteCompletionFile {
        input:
            outdir = base_out_dir + "/",
            keyfile =  MergePrimaryGenomeAlignedArrayElements.merged_bai
    }
}
