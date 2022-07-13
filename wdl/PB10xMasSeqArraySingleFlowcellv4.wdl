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

import "tasks/StringTie2.wdl"

import "tasks/TranscriptAnalysis/UMI_Tools.wdl" as UMI_TOOLS
import "tasks/TranscriptAnalysis/Postprocessing_Tasks.wdl" as TX_POST
import "tasks/TranscriptAnalysis/Preprocessing_Tasks.wdl" as TX_PRE

workflow PB10xMasSeqSingleFlowcellv4 {

    meta {
        description : "This workflow is designed to process data from the MASSeq v2 protocol and produce aligned reads that are ready for downstream analysis (e.g. transcript isoform identification).  It takes in a raw PacBio run folder location on GCS and produces a folder containing the aligned reads and other processed data."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String gcs_input_dir
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/PB10xMasSeqSingleFlowcellv4"

        File cell_barcode_whitelist = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/737K-august-2016.txt"

        # NOTE: Reference for un-split CCS reads:
        File ref_fasta =  "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
        File ref_fasta_index = "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai"
        File ref_fasta_dict = "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict"

        # NOTE: Reference for array elements:
        File transcriptome_ref_fasta =  "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.pc_transcripts.fa"
        File transcriptome_ref_fasta_index = "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.pc_transcripts.fa.fai"
        File transcriptome_ref_fasta_dict = "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.pc_transcripts.dict"

        File genome_annotation_gtf = "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.primary_assembly.annotation.gtf"

        File jupyter_template_static = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/MAS-seq_QC_report_template-static.ipynb"
        File workflow_dot_file = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/PB10xMasSeqArraySingleFlowcellv2.dot"

        File intervals_of_interest = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/gencode.37.TCR_intervals.tsv"
        String interval_overlap_name = "is_tcr_overlapping"

        File short_read_umis_tsv

        String starcode_extra_params = "--dist 2 --sphere"

        String expanded_cbc_tag = "CR"

        File? illumina_barcoded_bam

        # Default here is 0 because ccs uncorrected reads all seem to have RQ = -1.
        # All pathologically long reads also have RQ = -1.
        # This way we preserve the vast majority of the data, even if it has low quality.
        # We can filter it out at later steps.
        Float min_read_quality = 0.0
        Int max_reclamation_length = 60000

        Boolean is_SIRV_data = false
        String mas_seq_model = "mas_15_sc_10x5p_single_none"

        Int ccs_lev_dist = 2
        Int clr_lev_dist = 3

        Int ccs_umi_padding = 2
        Int clr_umi_padding = 2

        Int ccs_cbc_padding = 3
        Int clr_cbc_padding = 3

        # Add a suffix here for our out directory so we can label runs:
        String out_dir_suffix = ""

        String? sample_name
    }

    parameter_meta {
        gcs_input_dir : "Input folder on GCS in which to search for BAM files to process."
        gcs_out_root_dir : "Root output GCS folder in which to place results of this workflow."

        cell_barcode_whitelist : "Text file containing a whitelist of cell barcodes for the single-cell library prep."

        ref_fasta : "FASTA file containing the reference sequence to which the input data should be aligned before splitting into array elements."
        ref_fasta_index : "FASTA index file for the given ref_fasta file."
        ref_fasta_dict : "Sequence dictionary file for the given ref_fasta file."

        transcriptome_ref_fasta : "FASTA file containing the reference sequence to which the array elements should be aligned."
        transcriptome_ref_fasta_index : "FASTA index file for the given transcriptome_ref_fasta file."
        transcriptome_ref_fasta_dict : "Sequence dictionary file for the given transcriptome_ref_fasta file."

        genome_annotation_gtf : "Gencode GTF file containing genome annotations for the organism under study (usually humans).  This must match the given reference version and transcriptiome reference (usually hg38)."

        jupyter_template_static : "Jupyter notebook / ipynb file containing a template for the QC report which will contain static plots.  This should contain the same information as the jupyter_template_interactive file, but with static images."
        workflow_dot_file : "DOT file containing the representation of this WDL to be included in the QC reports.  This can be generated with womtool."

        intervals_of_interest : "[optional] An interval list file containing intervals to mark in the final anndata object as overlapping the transcripts.  Defaults to a T-cell receptor interval list."
        interval_overlap_name : "[optional] The name of the annotation to add to the final anndata object for the column containing the overlap flag for transcripts that overlap intervals in the given intervals_of_interest file.  Default: is_tcr_overlapping"

        illumina_barcoded_bam : "[optional] Illumina short reads file from a replicate of this same sample.  Used to perform cell barcode corrections."

        min_read_quality : "[optional] Minimum read quality for reads to have to be included in our data (Default: 0.0)."
        max_reclamation_length : "[optional] Maximum length (in bases) that a read can be to attempt to reclaim from CCS rejection (Default: 60000)."

        is_SIRV_data : "[optional] true if and only if the data in this sample are from the SIRV library prep.  false otherwise (Default: false)"
        mas_seq_model : "[optional] built-in mas-seq model to use (Default: mas15)"

        sample_name : "[optional] The name of the sample to associate with the data in this workflow."
    }

    # Create some runtime attributes that will force google to do network transfers really fast:
    RuntimeAttr fast_network_attrs = object {
        cpu_cores:  4,
        mem_gb:     32,
        disk_type:  "LOCAL",
        preemptible_tries:  0
    }
    RuntimeAttr fast_network_attrs_preemptible = object {
        cpu_cores:  4,
        mem_gb:     32,
        disk_type:  "LOCAL",
        preemptible_tries:  1
    }

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_01_WdlExecutionStartTimestamp { input: }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams as t_02_FindBams { input: gcs_input_dir = gcs_input_dir }
    call PB.FindZmwStatsJsonGz as t_03_FindZmwStatsJsonGz { input: gcs_input_dir = gcs_input_dir }

    # Check here if we found ccs bams or subread bams:
    Boolean use_subreads = t_02_FindBams.has_subreads
    Array[String] top_level_bam_files = if use_subreads then t_02_FindBams.subread_bams else t_02_FindBams.ccs_bams

    # Make sure we have **EXACTLY** one bam file to run on:
    if (length(top_level_bam_files) != 1) {
        call Utils.FailWithWarning as t_04_WARN1 { input: warning = "Error: Multiple BAM files found.  Cannot continue!" }
    }

    if (use_subreads) {
        call Utils.FailWithWarning as t_05_WARN2 { input: warning = "Error: This workflow now only supports data from the Sequel IIe." }
    }

    # Alias our bam file so we can work with it easier:
    File reads_bam = top_level_bam_files[0]

    call PB.GetPbReadGroupInfo as t_06_GetReadGroupInfo { input: gcs_bam_path = reads_bam }
    call PB.GetRunInfo as t_07_GetRunInfo { input: subread_bam = reads_bam }

    String SM  = select_first([sample_name, t_07_GetRunInfo.run_info["SM"]])
    String PL  = "PACBIO"
    String PU  = t_07_GetRunInfo.run_info["PU"]
    String DT  = t_07_GetRunInfo.run_info["DT"]
    String ID  = PU
    String DS  = t_07_GetRunInfo.run_info["DS"]
    String DIR = SM + "." + ID

    String RG_subreads  = "@RG\\tID:~{ID}.subreads\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
    String RG_consensus = "@RG\\tID:~{ID}.consensus\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
    String RG_array_elements = "@RG\\tID:~{ID}.array_elements\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

    # Check to see if we need to annotate our reads:
    call LONGBOW.CheckForAnnotatedArrayReads as t_08_CheckForAnnotatedReads {
        input:
            bam = reads_bam
    }

    File read_pbi = sub(reads_bam, ".bam$", ".bam.pbi")
    call PB.ShardLongReads as t_09_ShardLongReads {
        input:
            unaligned_bam = reads_bam,
            unaligned_pbi = read_pbi,
            prefix = SM + "_shard",
            num_shards = 300,
    }

    ## No more preemption on this sharding - takes too long otherwise.
    RuntimeAttr disable_preemption_runtime_attrs = object {
        preemptible_tries: 0
    }

    Array[String] tags_to_preserve =  [ "CB", "JB", "JC", "JD", "JF", "JX", "RC", "RG", "SG", "XA", "XB", "XC", "XF", "XM", "XN", "XQ", "XU", "YC", "YG", "YK", "YN", "YP", "YQ", "YS", "YV", "ZS", "ZU", "ec", "fn", "ic", "im", "is", "it", "np", "pz", "rn", "rq", "sn", "we", "ws", "zm" ]

    scatter (main_shard_index in range(length(t_09_ShardLongReads.unmapped_shards))) {
        File sharded_reads = t_09_ShardLongReads.unmapped_shards[main_shard_index]

        String fbmrq_prefix = basename(sharded_reads, ".bam")

        # Filter out the kinetics tags from PB files:
        call PB.RemoveKineticsTags as t_10_RemoveKineticsTags {
            input:
                bam = sharded_reads,
                prefix = SM + "_kinetics_removed"
        }

        # Handle setting up the things that we need for further processing of CCS-only reads:
        call PB.FindCCSReport as t_11_FindCCSReport {
            input:
                gcs_input_dir = gcs_input_dir
        }

        # 1 - filter the reads by the minimum read quality:
        call Utils.Bamtools as t_12_FilterS2EByMinReadQuality {
            input:
                bamfile = t_10_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_good_reads",
                cmd = "filter",
                args = '-tag "rq":">=' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # 1.5 - Get the "rejected" reads:
        call Utils.Bamtools as t_13_GetS2ECcsRejectedReads {
            input:
                bamfile = t_10_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_rejected_reads",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        #################################################################################################################
        # 2 - Get reads we can reclaim:
        call Utils.Bamtools as t_14_ExtractS2ECcsReclaimableReads {
            input:
                bamfile = t_10_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_reads_for_ccs_reclamation",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '" -length "<=' + max_reclamation_length + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        if ( ! t_08_CheckForAnnotatedReads.bam_has_annotations ) {
            # 3: Longbow annotate ccs reads
            call LONGBOW.Annotate as t_15_AnnotateS2ECCSReads {
                input:
                    reads = t_12_FilterS2EByMinReadQuality.bam_out,
                    model = mas_seq_model
            }
            # 4: Longbow annotate reclaimable reads
            call LONGBOW.Annotate as t_16_AnnotateS2EReclaimableReads {
                input:
                    reads = t_14_ExtractS2ECcsReclaimableReads.bam_out,
                    model = mas_seq_model
            }
        }

        File annotated_S2E_ccs_file = if t_08_CheckForAnnotatedReads.bam_has_annotations then t_12_FilterS2EByMinReadQuality.bam_out else select_first([t_15_AnnotateS2ECCSReads.annotated_bam])
        File annotated_S2E_reclaimable_file = if t_08_CheckForAnnotatedReads.bam_has_annotations then t_14_ExtractS2ECcsReclaimableReads.bam_out else select_first([t_16_AnnotateS2EReclaimableReads.annotated_bam])

        # 5: Longbow filter ccs annotated reads
        call LONGBOW.Filter as t_17_FilterS2ECCSReads {
            input:
                bam = annotated_S2E_ccs_file,
                prefix = SM + "_subshard",
                model = mas_seq_model
        }

        # 6: Longbow filter ccs reclaimable reads
        call LONGBOW.Filter as t_18_FilterS2EReclaimableReads {
            input:
                bam = annotated_S2E_reclaimable_file,
                prefix = SM + "_subshard",
                model = mas_seq_model
        }

        # 7: PBIndex CCS reads
        call PB.PBIndex as t_19_PbIndexS2ELongbowPassedCcsReads {
            input:
                bam = t_17_FilterS2ECCSReads.passed_reads
        }

        call PB.PBIndex as t_20_PbIndexS2ELongbowFailedCcsReads {
            input:
                bam = t_17_FilterS2ECCSReads.failed_reads
        }

        # 8: PBIndex reclaimable reads
        call PB.PBIndex as t_21_PbIndexS2ELongbowPassedReclaimedReads {
            input:
                bam = t_18_FilterS2EReclaimableReads.passed_reads
        }
        call PB.PBIndex as t_22_PbIndexS2ELongbowFailedReclaimableReads {
            input:
                bam = t_18_FilterS2EReclaimableReads.failed_reads
        }

        #####################################################################################################################
        # Now we have CCS and Reclaimed reads.
        #############################################

        # New pipeline steps:
        #     (1) minimap2 CCS reads in hifi mode
        #     (2) minimap2 CLR reads in shitty mode
        #     (3) merge 1 and 2 alignments
        #     <ADD FILTERING HERE>
        #     D2 sphere for cell barcode correction
        #     (4) annotate merged alignments with SQANTI3
        #     (5) cell barcode error correction with starcode (Levenshtein 1)
        #     (6) UMI error correction with umi-tools (Levenshtein 2), grouping by SQANTI3
        #     (7) generate count matrix

        # Shard our CCS reads into smaller problems to do work on array elements:
        call PB.ShardLongReads as t_23_ShardS2ECcsLongbowPassedReads {
            input:
                unaligned_bam = t_17_FilterS2ECCSReads.passed_reads,
                unaligned_pbi = t_19_PbIndexS2ELongbowPassedCcsReads.pbindex,
                prefix = SM + "_ccs_reads_shard_" + main_shard_index,
                num_shards = 10,
        }

        scatter (ccsi1 in range(length(t_23_ShardS2ECcsLongbowPassedReads.unmapped_shards))) {
            File s2e_ccs_longbow_passed_shard = t_23_ShardS2ECcsLongbowPassedReads.unmapped_shards[ccsi1]

            # Segment CCS reads into array elements:
            call LONGBOW.Segment as t_24_SegmentS2ECcsReads {
                input:
                    annotated_reads = s2e_ccs_longbow_passed_shard,
                    prefix = SM + "_ccs_array_elements_shard_" + main_shard_index + "_subshard_" + ccsi1,
                    extra_args = "-b",
                    model = mas_seq_model
            }

            # Now remove all -END reads:
            # We do this here to speed up all other calculations.
            call Utils.RemoveMasSeqTruncatedReads as t_25_RemoveMasSeqTruncatedReadsFromCcsReads {
                input:
                    bam_file = t_24_SegmentS2ECcsReads.segmented_bam
            }

            # Now that we've annotated the reads, we can pad the UMIs by a couple of bases to aid in the deduping:
            call LONGBOW.Pad as t_26_LongbowPadCCSArrayElementUMIs {
                input:
                    reads = t_25_RemoveMasSeqTruncatedReadsFromCcsReads.bam,
                    model = mas_seq_model,
                    tag_to_expand = "ZU",
                    padding = ccs_umi_padding,
                    prefix = SM + "_ccs_array_elements_annotated_umi_padded_shard_" + main_shard_index + "_subshard_" + ccsi1
            }

            call LONGBOW.Pad as t_27_LongbowPadCCSArrayElementCBCs {
                input:
                    reads = t_26_LongbowPadCCSArrayElementUMIs.padded_tag_bam,
                    model = mas_seq_model,
                    tag_to_expand = "CR",
                    new_tag_dest = expanded_cbc_tag,
                    padding = ccs_cbc_padding,
                    prefix = SM + "_ccs_array_elements_annotated_cbc_padded_shard_" + main_shard_index + "_subshard_" + ccsi1
            }

            # Now we should correct our barcodes based on the whitelist:
            call LONGBOW.Correct as t_28_LongbowCorrectCCSCorrectedArrayElementCBCs {
                input:
                    reads = t_27_LongbowPadCCSArrayElementCBCs.padded_tag_bam,
                    barcode_allow_list = cell_barcode_whitelist,
                    model = mas_seq_model,
                    ccs_lev_dist_threshold = ccs_lev_dist,
                    clr_lev_dist_threshold = clr_lev_dist,
                    prefix = SM + "_ccs_array_elements_annotated_padded_cbc_corrected_shard_" + main_shard_index + "_subshard_" + ccsi1,
                    raw_barcode_tag = expanded_cbc_tag,
                    corrected_barcode_tag = "CB",
            }

            call TENX.AdjustUmiSequenceWithAdapterAlignment as t_29_AdjustCCSUMIs {
                input:
                    bam = t_28_LongbowCorrectCCSCorrectedArrayElementCBCs.corrected_barcodes_bam,
                    short_read_umis = short_read_umis_tsv,
                    prefix = SM + "_ccs_array_elements_annotated_padded_cbc_corrected_UMI_adjusted_shard_" + main_shard_index + "_subshard_" + ccsi1
            }

            call LONGBOW.Extract as t_30_LongbowExtractCcsArrayElements {
                input:
                    bam = t_29_AdjustCCSUMIs.output_bam,
                    prefix = SM + "_ccs_array_elements_cbc_umi_padded_extracted_shard_" + main_shard_index + "_subshard_" + ccsi1
            }
        }

        # Merge Filtered CCS reads together:
        call Utils.MergeBams as t_31_MergeCCSArrayElementShards {
            input:
                bams = t_24_SegmentS2ECcsReads.segmented_bam,
                prefix = SM + "_ccs_array_elements_shard_" + main_shard_index
        }

        # Merge CCS Barcode Conf files:
        call Utils.MergeFiles as t_32_MergeCCSBarcodeConfShards {
            input:
                files_to_merge = t_24_SegmentS2ECcsReads.barcode_conf_file,
                merged_file_name = SM + "_ccs_array_element_barcode_confs_shard_" + main_shard_index + ".txt"
        }

        # Merge Filtered CCS reads with no ends together:
        call Utils.MergeBams as t_33_MergeCCSArrayElementsNonTruncatedShards {
            input:
                bams = t_25_RemoveMasSeqTruncatedReadsFromCcsReads.bam,
                prefix = SM + "_ccs_array_elements_no_ends_shard_" + main_shard_index
        }

        # Merge UMI-Padded CCS Array Elements:
        call Utils.MergeBams as t_34_MergeCCSArrayElementsUmiPaddedShards {
            input:
                bams = t_26_LongbowPadCCSArrayElementUMIs.padded_tag_bam,
                prefix = SM + "_ccs_array_elements_no_ends_umi_padded_shard_" + main_shard_index
        }

        # Merge CBC-UMI-Padded CCS Array Elements:
        call Utils.MergeBams as t_35_MergeCCSArrayElementsUmiCbcPaddedShards {
            input:
                bams = t_27_LongbowPadCCSArrayElementCBCs.padded_tag_bam,
                prefix = SM + "_ccs_array_elements_no_ends_cbc_umi_padded_shard_" + main_shard_index
        }

        # Merge Corrected CBC CCS Array Elements:
        call Utils.MergeBams as t_36_MergeCCSArrayElementsUmiCbcPaddedCbcCorrectedShards {
            input:
                bams = t_28_LongbowCorrectCCSCorrectedArrayElementCBCs.corrected_barcodes_bam,
                prefix = SM + "_ccs_array_elements_no_ends_cbc_umi_padded_shard_" + main_shard_index
        }

        call Utils.MergeBams as t_37_MergeLongbowPaddedCBCUncorrectableCCSArrayElementsShards {
            input:
                bams = t_28_LongbowCorrectCCSCorrectedArrayElementCBCs.uncorrected_barcodes_bam,
                prefix = SM + "_ccs_array_elements_aligned_annotated_padded_CBC_uncorrectable_shard_" + main_shard_index
        }

        call Utils.MergeBams as t_38_MergeLongbowPaddedCBCCorrectedCCSArrayElementsShards {
            input:
                bams = t_29_AdjustCCSUMIs.output_bam,
                prefix = SM + "_ccs_array_elements_aligned_annotated_padded_CBC_corrected_shard_" + main_shard_index
        }

        call Utils.MergeBams as t_39_MergeLongbowExtractedCcsArrayElements {
            input:
                bams = t_30_LongbowExtractCcsArrayElements.extracted_bam,
                prefix = SM + "_ccs_array_elements_cbc_umi_padded_extracted_shard_" + main_shard_index
        }

        #####################
        # CCS Reclaimed / CLR:

        # Shard our CCS reclaimed reads into smaller problems to do work on array elements:
        call PB.ShardLongReads as t_40_ShardS2ECcsReclaimedReads {
            input:
                unaligned_bam = t_18_FilterS2EReclaimableReads.passed_reads,
                unaligned_pbi = t_21_PbIndexS2ELongbowPassedReclaimedReads.pbindex,
                prefix = SM + "_ccs_reclaimed_reads_shard_" + main_shard_index,
                num_shards = 10,
        }

        scatter (cri1 in range(length(t_40_ShardS2ECcsReclaimedReads.unmapped_shards))) {
            File s2e_ccs_reclaimed_shard = t_40_ShardS2ECcsReclaimedReads.unmapped_shards[cri1]

            # Segment Reclaimed reads into array elements:
            call LONGBOW.Segment as t_41_SegmentS2ECcsReclaimedReads {
                input:
                    annotated_reads = s2e_ccs_reclaimed_shard,
                    prefix = SM + "_ccs_reclaimed_array_elements_shard_" + main_shard_index + "_subshard_" + cri1,
                    extra_args = "-b",
                    model = mas_seq_model
            }

            # Now remove all -END reads:
            # We do this here to speed up all other calculations.
            call Utils.RemoveMasSeqTruncatedReads as t_42_RemoveMasSeqTruncatedReadsFromCcsReclaimedReads {
                input:
                    bam_file = t_41_SegmentS2ECcsReclaimedReads.segmented_bam
            }

            # Now that we've annotated the reads, we can pad the UMIs by a couple of bases to aid in the deduping:
            call LONGBOW.Pad as t_43_LongbowPadCcsReclaimedArrayElementUMIs {
                input:
                    reads = t_42_RemoveMasSeqTruncatedReadsFromCcsReclaimedReads.bam,
                    model = mas_seq_model,
                    tag_to_expand = "ZU",
                    padding = ccs_umi_padding,
                    prefix = SM + "_ccs_reclaimed_array_elements_annotated_umi_padded_shard_" + main_shard_index + "_subshard_" + cri1
            }

            call LONGBOW.Pad as t_44_LongbowPadCcsReclaimedArrayElementCBCs {
                input:
                    reads = t_43_LongbowPadCcsReclaimedArrayElementUMIs.padded_tag_bam,
                    model = mas_seq_model,
                    tag_to_expand = "CR",
                    new_tag_dest = expanded_cbc_tag,
                    padding = ccs_cbc_padding,
                    prefix = SM + "_ccs_reclaimed_array_elements_annotated_cbc_padded_shard_" + main_shard_index + "_subshard_" + cri1
            }

            # Now we should correct our barcodes based on the whitelist:
            call LONGBOW.Correct as t_45_LongbowCorrectCcsReclaimedArrayElementCBCs {
                input:
                    reads = t_44_LongbowPadCcsReclaimedArrayElementCBCs.padded_tag_bam,
                    barcode_allow_list = cell_barcode_whitelist,
                    model = mas_seq_model,
                    ccs_lev_dist_threshold = ccs_lev_dist,
                    clr_lev_dist_threshold = clr_lev_dist,
                    prefix = SM + "_ccs_reclaimed_array_elements_annotated_padded_cbc_corrected_shard_" + main_shard_index + "_subshard_" + cri1,
                    raw_barcode_tag = expanded_cbc_tag,
                    corrected_barcode_tag = "CB",
            }

            call TENX.AdjustUmiSequenceWithAdapterAlignment as t_46_AdjustCcsReclaimedUMIs {
                input:
                    bam = t_45_LongbowCorrectCcsReclaimedArrayElementCBCs.corrected_barcodes_bam,
                    short_read_umis = short_read_umis_tsv,
                    prefix = SM + "_ccs_reclaimed_array_elements_annotated_padded_cbc_corrected_UMI_adjusted_shard_" + main_shard_index + "_subshard_" + cri1
            }

            call LONGBOW.Extract as t_47_LongbowExtractCcsReclaimedArrayElements {
                input:
                    bam = t_46_AdjustCcsReclaimedUMIs.output_bam,
                    prefix = SM + "_ccs_reclaimed_array_elements_cbc_umi_padded_extracted_shard_" + main_shard_index + "_subshard_" + cri1
            }
        }

        # Merge Filtered CCS Reclaimed reads together:
        call Utils.MergeBams as t_48_MergeCCSReclaimedArrayElementShards {
            input:
                bams = t_41_SegmentS2ECcsReclaimedReads.segmented_bam,
                prefix = SM + "_ccs_reclaimed_array_elements_shard_" + main_shard_index
        }

        # Merge CCS Barcode Conf files:
        call Utils.MergeFiles as t_49_MergeCCSReclaimedBarcodeConfShards {
            input:
                files_to_merge = t_41_SegmentS2ECcsReclaimedReads.barcode_conf_file,
                merged_file_name = SM + "_ccs_reclaimed_array_element_barcode_confs_shard_" + main_shard_index + ".txt"
        }

        # Merge Filtered CCS reads with no ends together:
        call Utils.MergeBams as t_50_MergeCCSReclaimedArrayElementsNonTruncatedShards {
            input:
                bams = t_42_RemoveMasSeqTruncatedReadsFromCcsReclaimedReads.bam,
                prefix = SM + "_ccs_reclaimed_array_elements_no_ends_shard_" + main_shard_index
        }

        # Merge UMI-Padded CCS Array Elements:
        call Utils.MergeBams as t_51_MergeCCSReclaimedArrayElementsUmiPaddedShards {
            input:
                bams = t_43_LongbowPadCcsReclaimedArrayElementUMIs.padded_tag_bam,
                prefix = SM + "_ccs_reclaimed_array_elements_no_ends_umi_padded_shard_" + main_shard_index
        }

        # Merge CBC-UMI-Padded CCS Array Elements:
        call Utils.MergeBams as t_52_MergeCCSReclaimedArrayElementsUmiCbcPaddedShards {
            input:
                bams = t_44_LongbowPadCcsReclaimedArrayElementCBCs.padded_tag_bam,
                prefix = SM + "_ccs_reclaimed_array_elements_no_ends_cbc_umi_padded_shard_" + main_shard_index
        }

        # Merge Corrected CBC CCS Array Elements:
        call Utils.MergeBams as t_53_MergeCCSReclaimedArrayElementsUmiCbcPaddedCbcCorrectedShards {
            input:
                bams = t_45_LongbowCorrectCcsReclaimedArrayElementCBCs.corrected_barcodes_bam,
                prefix = SM + "_ccs_reclaimed_array_elements_no_ends_cbc_umi_padded_shard_" + main_shard_index
        }

        call Utils.MergeBams as t_54_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElementsShards {
            input:
                bams = t_45_LongbowCorrectCcsReclaimedArrayElementCBCs.uncorrected_barcodes_bam,
                prefix = SM + "_ccs_reclaimed_array_elements_aligned_annotated_padded_CBC_uncorrectable_shard_" + main_shard_index
        }

        call Utils.MergeBams as t_55_MergeLongbowPaddedCBCCorrectedCCSReclaimedArrayElementsShards {
            input:
                bams = t_46_AdjustCcsReclaimedUMIs.output_bam,
                prefix = SM + "_ccs_reclaimed_array_elements_aligned_annotated_padded_CBC_corrected_shard_" + main_shard_index
        }

        call Utils.MergeBams as t_56_MergeLongbowExtractedCcsReclaimedArrayElements {
            input:
                bams = t_47_LongbowExtractCcsReclaimedArrayElements.extracted_bam,
                prefix = SM + "_ccs_reclaimed_array_elements_cbc_umi_padded_extracted_shard_" + main_shard_index
        }

        ###############

        # Now align the array elements with their respective alignment presets.
        # NOTE: We use the non-truncated reads because we only want the good stuff.

        # Align CCS reads to the genome:
        call AR.Minimap2 as t_57_AlignCCSArrayElementsToGenome {
            input:
                reads      = [ t_39_MergeLongbowExtractedCcsArrayElements.merged_bam ],
                ref_fasta  = ref_fasta,
                tags_to_preserve = tags_to_preserve,
                map_preset = "splice:hq",
                prefix = SM + "_ccs_array_elements_extracted_aligned_shard_"+ main_shard_index,
                runtime_attr_override = object { mem_gb: 32 }
        }

        # Align Reclaimed reads to the genome:
        call AR.Minimap2 as t_58_AlignReclaimedArrayElementsToGenome {
            input:
                reads      = [ t_56_MergeLongbowExtractedCcsReclaimedArrayElements.merged_bam ],
                ref_fasta  = ref_fasta,
                tags_to_preserve = tags_to_preserve,
                map_preset = "splice",
                prefix = SM + "_ccs_reclaimed_array_elements_extracted_aligned_shard_" + main_shard_index,
                runtime_attr_override = object { mem_gb: 32 }
        }
    }

    # Merge the arrays:
    call Utils.MergeBams as t_59_MergeCCSLongbowAnnotatedArrayReads {
        input:
            bams = annotated_S2E_ccs_file,
            prefix = SM + "_ccs_array_reads_longbow_annotated"
    }
    call PB.PBIndex as t_60_PbIndexMergedCCSLongbowAnnotatedArrayReads { input: bam = t_59_MergeCCSLongbowAnnotatedArrayReads.merged_bam }

    call Utils.MergeBams as t_61_MergeCCSReclaimableLongbowAnnotatedArrayReads {
        input:
            bams = annotated_S2E_reclaimable_file,
            prefix = SM + "_ccs_reclaimable_array_reads_longbow_annotated"
    }
    call PB.PBIndex as t_62_PbIndexMergedCCSReclaimableLongbowAnnotatedArrayReads { input: bam = t_61_MergeCCSReclaimableLongbowAnnotatedArrayReads.merged_bam }

    call Utils.MergeBams as t_63_MergeCCSLongbowPassedArrayReads {
        input:
            bams = t_17_FilterS2ECCSReads.passed_reads,
            prefix = SM + "_ccs_array_reads_longbow_passed"
    }
    call PB.PBIndex as t_64_PbIndexMergedCCSLongbowPassedArrayReads { input: bam = t_63_MergeCCSLongbowPassedArrayReads.merged_bam }

    call Utils.MergeBams as t_65_MergeCCSLongbowFailedArrayReads {
        input:
            bams = t_17_FilterS2ECCSReads.failed_reads,
            prefix = SM + "_ccs_array_reads_longbow_failed"
    }
    call PB.PBIndex as t_66_PbIndexMergedCCSLongbowFailedArrayReads { input: bam = t_65_MergeCCSLongbowFailedArrayReads.merged_bam }

    call Utils.MergeBams as t_67_MergeCCSReclaimedArrayReads {
        input:
            bams = t_18_FilterS2EReclaimableReads.passed_reads,
            prefix = SM + "_ccs_reclaimed_array_reads_longbow_passed"
    }
    call PB.PBIndex as t_68_PbIndexMergedCCSReclaimedArrayReads { input: bam = t_67_MergeCCSReclaimedArrayReads.merged_bam }

    call Utils.MergeBams as t_69_MergeCCSUnreclaimableArrayReads {
        input:
            bams = t_18_FilterS2EReclaimableReads.failed_reads,
            prefix = SM + "_ccs_unreclaimable_array_reads_longbow_failed"
    }
    call PB.PBIndex as t_70_PbIndexMergedCCSUnreclaimableReclaimedArrayReads { input: bam = t_69_MergeCCSUnreclaimableArrayReads.merged_bam }

    call Utils.MergeBams as t_71_MergeLongbowPassedReads {
        input:
            bams = flatten([t_17_FilterS2ECCSReads.passed_reads, t_18_FilterS2EReclaimableReads.passed_reads]),
            prefix = SM + "_longbow_passed_array_reads"
    }
    call PB.PBIndex as t_72_PbIndexMergedLongbowPassingReads { input: bam = t_71_MergeLongbowPassedReads.merged_bam }

    call Utils.MergeBams as t_73_MergeLongbowFailedReads {
        input:
            bams = flatten([t_17_FilterS2ECCSReads.failed_reads, t_18_FilterS2EReclaimableReads.failed_reads]),
            prefix = SM + "_longbow_failed_array_reads"
    }
    call PB.PBIndex as t_74_PbIndexMergedLongbowFailedReads { input: bam = t_73_MergeLongbowFailedReads.merged_bam }

    call Utils.MergeBams as t_75_MergeAllLongbowAnnotatedReads {
        input:
            bams = flatten([annotated_S2E_ccs_file, annotated_S2E_reclaimable_file]),
            prefix = SM + "_all_longbow_annotated_array_reads"
    }
    call PB.PBIndex as t_76_PbIndexMergedAllLongbowAnnotatedReads { input: bam = t_75_MergeAllLongbowAnnotatedReads.merged_bam }

    # Merge Filtered CCS array elements together:
    call Utils.MergeBams as t_77_MergeCCSArrayElements {
        input:
            bams = t_31_MergeCCSArrayElementShards.merged_bam,
            prefix = SM + "_ccs_array_elements"
    }

    # Merge Filtered CCS Reclaimed array elements together:
    call Utils.MergeBams as t_78_MergeCCSReclaimedArrayElements {
        input:
            bams = t_48_MergeCCSReclaimedArrayElementShards.merged_bam,
            prefix = SM + "_ccs_reclaimed_array_elements"
    }

    # Merge Filtered (and end removed) CCS array elements together:
    call Utils.MergeBams as t_79_MergeCCSArrayElementsNonTruncated {
        input:
            bams = t_33_MergeCCSArrayElementsNonTruncatedShards.merged_bam,
            prefix = SM + "_ccs_array_elements_non_truncated"
    }

    # Merge Filtered (and end removed) CCS Reclaimed array elements together:
    call Utils.MergeBams as t_80_MergeCCSReclaimedArrayElementsNonTruncated {
        input:
            bams = t_50_MergeCCSReclaimedArrayElementsNonTruncatedShards.merged_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_non_truncated"
    }

    call Utils.MergeBams as t_81_MergeAllArrayElementsNonTruncated {
        input:
            bams = flatten([t_33_MergeCCSArrayElementsNonTruncatedShards.merged_bam, t_50_MergeCCSReclaimedArrayElementsNonTruncatedShards.merged_bam]),
            prefix = SM + "_array_elements_non_truncated"
    }

    # CCS Corrected:

    # Merge UMI-Padded CCS Array Elements:
    call Utils.MergeBams as t_82_MergeCCSArrayElementsUmiPadded {
        input:
            bams = t_34_MergeCCSArrayElementsUmiPaddedShards.merged_bam,
            prefix = SM + "_ccs_array_elements_no_ends_umi_padded"
    }

    # Merge CBC-UMI-Padded CCS Array Elements:
    call Utils.MergeBams as t_83_MergeCCSArrayElementsUmiCbcPadded {
        input:
            bams = t_35_MergeCCSArrayElementsUmiCbcPaddedShards.merged_bam,
            prefix = SM + "_ccs_array_elements_no_ends_cbc_umi_padded"
    }

    # Merge Corrected CBC CCS Array Elements:
    call Utils.MergeBams as t_84_MergeCCSArrayElementsUmiCbcPaddedCbcCorrected {
        input:
            bams = t_36_MergeCCSArrayElementsUmiCbcPaddedCbcCorrectedShards.merged_bam,
            prefix = SM + "_ccs_array_elements_no_ends_cbc_umi_padded"
    }

    # Merge Corrected CBC CCS Array Elements:
    call Utils.MergeBams as t_85_MergeLongbowPaddedCBCCorrectedCCSArrayElements {
        input:
            bams = t_38_MergeLongbowPaddedCBCCorrectedCCSArrayElementsShards.merged_bam,
            prefix = SM + "_ccs_array_elements_aligned_annotated_padded_CBC_corrected"
    }

    # Merge uncorrectable CBC CCS Array Elements:
    call Utils.MergeBams as t_86_MergeLongbowPaddedCBCUncorrectableCCSArrayElements {
        input:
            bams = t_37_MergeLongbowPaddedCBCUncorrectableCCSArrayElementsShards.merged_bam,
            prefix = SM + "_ccs_array_elements_aligned_annotated_padded_CBC_uncorrectable"
    }

    call Utils.MergeBams as t_87_MergeLongbowPaddedCBCCorrectedExtractedCCSArrayElements {
        input:
            bams = t_39_MergeLongbowExtractedCcsArrayElements.merged_bam,
            prefix = SM + "_ccs_array_elements_cbc_umi_padded_extracted"
    }

    # RECLAIMED:

    # Merge UMI-Padded CCS Reclaimed Array Elements:
    call Utils.MergeBams as t_88_MergeCCSReclaimedArrayElementsUmiPadded {
        input:
            bams = t_51_MergeCCSReclaimedArrayElementsUmiPaddedShards.merged_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_no_ends_umi_padded"
    }

    # Merge CBC-UMI-Padded CCS Reclaimed Array Elements:
    call Utils.MergeBams as t_89_MergeCCSReclaimedArrayElementsUmiCbcPadded {
        input:
            bams = t_52_MergeCCSReclaimedArrayElementsUmiCbcPaddedShards.merged_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_no_ends_cbc_umi_padded"
    }

    # Merge Corrected CBC CCS Reclaimed Array Elements:
    call Utils.MergeBams as t_90_MergeCCSReclaimedArrayElementsUmiCbcPaddedCbcCorrected {
        input:
            bams = t_53_MergeCCSReclaimedArrayElementsUmiCbcPaddedCbcCorrectedShards.merged_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_no_ends_cbc_umi_padded"
    }

    # Merge uncorrectable CBC CCS Reclaimed Array Elements:
    call Utils.MergeBams as t_91_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElements {
        input:
            bams = t_54_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElementsShards.merged_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_aligned_annotated_padded_CBC_uncorrectable"
    }

    call Utils.MergeBams as t_92_MergeLongbowPaddedCBCCorrectedExtractedCCSReclaimedArrayElements {
        input:
            bams = t_56_MergeLongbowExtractedCcsReclaimedArrayElements.merged_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_cbc_umi_padded_extracted"
    }

    ##############################################################################################################


    # Merge Aligned CCS array elements together:
    call Utils.MergeBams as t_93_MergeAlignedCCSArrayElements {
        input:
            bams = t_57_AlignCCSArrayElementsToGenome.aligned_bam,
            prefix = SM + "_ccs_array_elements_padded_aligned"
    }

    # Merge Aligned CCS Reclaimed array elements together:
    call Utils.MergeBams as t_94_MergeAlignedCCSReclaimedArrayElements {
        input:
            bams = t_58_AlignReclaimedArrayElementsToGenome.aligned_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_padded_aligned"
    }

    call Utils.MergeBams as t_95_MergeAllAlignedArrayElementsNonTruncated {
        input:
            bams = flatten([t_57_AlignCCSArrayElementsToGenome.aligned_bam, t_58_AlignReclaimedArrayElementsToGenome.aligned_bam]),
            prefix = SM + "_array_elements_padded_aligned"
    }

    # Merge CCS Barcode Conf files:
    call Utils.MergeFiles as t_96_MergeAllCCSBarcodeConfShards {
        input:
            files_to_merge = t_32_MergeCCSBarcodeConfShards.merged_file,
            merged_file_name = SM + "_ccs_array_element_barcode_confs.txt"
    }

    # Merge CCS Barcode Conf files:
    call Utils.MergeFiles as t_97_MergeAllCCSReclaimedBarcodeConfShards {
        input:
            files_to_merge = t_49_MergeCCSReclaimedBarcodeConfShards.merged_file,
            merged_file_name = SM + "_ccs_reclaimed_array_element_barcode_confs.txt"
    }

    ##########################################################################################################################
    #############################################################
    #
    # FILTER THE READS HERE!!!!!!!!!!!!!!!!111!1!11
    #
    ############################################################
    ##########################################################################################################################

#    	• Post-Alignment filters:
#		○ Remove Flags:
#			§ MQ0
#			§ Supplementary
#			§ Secondary
#			§ Unmapped
#		○ Read length
#			§ Keep: <15Kb
#		○ Spread on the reference
#			§ Some fraciton of ref length
#			§ Or max N (splice gap)
#			§ Overlapping multiple genes - remove it.
#				□ Use funcotate segments or similar
#				□ Whitelist genes with overlapping exons
#				□ Bedtools intersection
#					® THOUGH, occasionally it gives different results
#					® Might not be a problem.
#		○ End soft/hard Clipping
#			§ L or R end of read
#               1000?

    RuntimeAttr filterReadsAttrs = object {
        preemptible_tries: 0
    }

    # Remove unmapped, secondary, supplementary, mq0, length > 15kb, end clips > 1kb
    call Utils.FilterMasSeqReadsWithGatk as t_98_AlignmentFilterForCcsArrayElements {
        input:
            bam_file = t_93_MergeAlignedCCSArrayElements.merged_bam,
            bam_index = t_93_MergeAlignedCCSArrayElements.merged_bai,
            prefix = SM + "_CCS_ArrayElements_Annotated_Aligned_PrimaryOnly",
            runtime_attr_override = filterReadsAttrs
    }

    call Utils.FilterMasSeqReadsWithGatk as t_99_AlignmentFilterForReclaimedArrayElements {
        input:
            bam_file = t_94_MergeAlignedCCSReclaimedArrayElements.merged_bam,
            bam_index = t_94_MergeAlignedCCSReclaimedArrayElements.merged_bai,
            prefix = SM + "_Reclaimed_ArrayElements_Annotated_Aligned_PrimaryOnly",
            runtime_attr_override = filterReadsAttrs
    }

    #########################################################################################################################
    ##########################################################################################################################

    # Mehrtash's suggestion for paper
    # 		• run stringtie2
    #		• filter stringtie2 annotations, get rid of super low TPM annotations
    #		polish stringtie2 annotations (e.g. if a transcript is an extension of a GENCODE transcript, propagate the name, if a “novel” stringtie2 gene overlaps with a previously annotated gene > 95%, propagate the name; otherwise, ignore)
    #		run TALON w/ polished stringtie2 annotations
    #       ignore NIC and NNC TALON transcripts (which should be VERY few), only focus on exact matches and ISM

    # Mehrtash's latest suggestion for paper:
    # Running StringTie2 on our sample w/ GENCODE as reference
    # Running TALON on GENCODE "reads" (we need to cook up a bam file from the GENCODE gtf) with a database initialized with the StringTie2 gtf
    # Running TALON on our sample with a darabase initialized with the StringTie2 gtf

    # Merge all alignments together:
    call Utils.MergeBams as t_100_MergeAllAlignedAndFilteredArrayElements {
        input:
            bams = [t_98_AlignmentFilterForCcsArrayElements.bam, t_99_AlignmentFilterForReclaimedArrayElements.bam],
            prefix = SM + "_all_array_elements_aligned_for_txome_discovery"
    }

    call StringTie2.Quantify as t_101_ST2_Quant {
        input:
            aligned_bam = t_100_MergeAllAlignedAndFilteredArrayElements.merged_bam,
            aligned_bai = t_100_MergeAllAlignedAndFilteredArrayElements.merged_bai,
            gtf = genome_annotation_gtf,
            keep_retained_introns = false,
            prefix = SM + "_StringTie2_Quantify",
    }

    call StringTie2.ExtractTranscriptSequences as t_102_ST2_ExtractTranscriptSequences  {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_index,
            gtf = t_101_ST2_Quant.st_gtf,
            prefix = SM + "_StringTie2_ExtractTranscriptSequences",
    }

    call StringTie2.CompareTranscriptomes as t_103_ST2_CompareTranscriptomes {
        input:
            guide_gtf = genome_annotation_gtf,
            new_gtf = t_101_ST2_Quant.st_gtf,
            prefix = SM + "_StringTie2_CompareTranscriptome",
    }

    ##########################################################################################################################
    ##########################################################################################################################

    #################
    # Here we restore the original read names to the bam because we're hashing them with Longbow.segment:

    # Restore original read names to CCS reads:
    call TX_PRE.RestoreOriginalReadNames as t_104_RestoreCcsOriginalReadNames {
        input:
            bam = t_98_AlignmentFilterForCcsArrayElements.bam,
            prefix =  SM + "_CCS_cbc_annotated_array_elements_padded_original_names"
    }

    # Restore original read names to CLR reads:
    call TX_PRE.RestoreOriginalReadNames as t_105_RestoreClrOriginalReadNames {
        input:
            bam = t_99_AlignmentFilterForReclaimedArrayElements.bam,
            prefix =  SM + "_CLR_cbc_annotated_array_elements_padded_original_names"
    }

    # Merge Aligned CCS and Reclaimed reads together:
    call Utils.MergeBams as t_106_MergeAllAnnotatedArrayElementsWithOriginalNames {
        input:
            bams = [t_104_RestoreCcsOriginalReadNames.bam_out, t_105_RestoreClrOriginalReadNames.bam_out],
            prefix = SM + "_all_cbc_annotated_array_elements_padded_original_names"
    }

    #################
    # Now we have to split the reads again, process them into gff files, run gffcompare and then aggregate the results in a graph

    # We can actually compare the references without needing to scatter:
    call TX_PRE.GffCompare as t_107_GffCompareStringtie2toGencode {
        input:
            gff_ref = t_101_ST2_Quant.st_gtf,
            gff_query = genome_annotation_gtf,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
    }
    call TX_PRE.GffCompare as t_108_GffCompareGencodetoStringtie2 {
        input:
            gff_ref = genome_annotation_gtf,
            gff_query = t_101_ST2_Quant.st_gtf,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
    }

    # Split by contig:
    call TX_PRE.SplitBamByContig as t_109_SplitArrayElementsByContig {
        input:
            bam = t_106_MergeAllAnnotatedArrayElementsWithOriginalNames.merged_bam,
            prefix = SM + "_all_cbc_annotated_array_elements_padded_original_names"
    }

    # For each contig:
    scatter (i in range(length(t_109_SplitArrayElementsByContig.contig_bams))) {

        File contig_bam = t_109_SplitArrayElementsByContig.contig_bams[i]
        String contig_name = t_109_SplitArrayElementsByContig.contig_names[i]

        # Create a GFF file:
        call TX_PRE.ConvertSplicedBamToGff as t_110_ConvertSplicedBamToGff {
            input:
                bam = contig_bam
        }

        # Compare GFF files:
        call TX_PRE.GffCompare as t_111_GffCompareStringtie2toMasSeqReads {
            input:
                gff_ref = t_101_ST2_Quant.st_gtf,
                gff_query = t_110_ConvertSplicedBamToGff.gff,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
        }

        call TX_PRE.GffCompare as t_112_GffCompareGencodetoMasSeqReads {
            input:
                gff_ref = genome_annotation_gtf,
                gff_query = t_110_ConvertSplicedBamToGff.gff,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
        }

        # Create the comparison graph and tsv files:
        call TX_POST.QuantifyGffComparison as t_113_QuantifyGffComparison {
            input:
                genome_gtf = genome_annotation_gtf,
                st2_gencode_refmap = t_107_GffCompareStringtie2toGencode.refmap,
                st2_gencode_tmap = t_107_GffCompareStringtie2toGencode.tmap,
                st2_read_refmap = t_111_GffCompareStringtie2toMasSeqReads.refmap,
                st2_read_tmap = t_111_GffCompareStringtie2toMasSeqReads.tmap,
                gencode_st2_refmap = t_108_GffCompareGencodetoStringtie2.refmap,
                gencode_st2_tmap = t_108_GffCompareGencodetoStringtie2.tmap,
                gencode_read_refmap = t_112_GffCompareGencodetoMasSeqReads.refmap,
                gencode_read_tmap = t_112_GffCompareGencodetoMasSeqReads.tmap,
                prefix = SM + "_all_cbc_annotated_array_elements_padded_" + contig_name
        }
    }

    # Merge our tx equivalance classes assignments and eq classes:
    call TX_POST.CombineEqClassFiles as t_114_CombineEqClassFiles {
        input:
            gene_eq_class_definitions = t_113_QuantifyGffComparison.gene_eq_class_labels_file,
            gene_assignment_files = t_113_QuantifyGffComparison.gene_assignments_file,
            equivalence_class_definitions = t_113_QuantifyGffComparison.tx_equivalence_class_labels_file,
            equivalence_classes = t_113_QuantifyGffComparison.tx_equivalence_class_file,
            prefix = SM + "_all_cbc_annotated_array_elements_padded"
    }

    ############################################################
    # Quantify Transcripts:
    ##########

    # Use old quant method here as a baseline for comparison:
    call TX_POST.CopyEqClassInfoToTag as t_115_CopyEqClassInfoToTag {
        input:
            bam = t_106_MergeAllAnnotatedArrayElementsWithOriginalNames.merged_bam,
            eq_class_file = t_114_CombineEqClassFiles.combined_tx_eq_class_assignments,
            prefix = SM + "_annotated_array_elements_for_quant_with_gene_names"
    }

    call TX_PRE.CorrectUmisWithSetCover as t_116_CorrectUmisWithSetCover {
        input:
            bam = t_115_CopyEqClassInfoToTag.bam_out,
            prefix = SM + "_annotated_array_elements_for_quant_with_gene_names"
    }

    # Because of how we're doing things, we need to pull out the CCS and CCS Reclaimed reads from the output of the
    # set cover correction:
    call Utils.Bamtools as t_117_GetCcsCorrectedReadsWithCorrectedUmis {
        input:
            bamfile = t_116_CorrectUmisWithSetCover.corrected_umi_reads,
            prefix = SM + "_annotated_array_elements_for_quant_with_gene_names.corrected_umis.CCS",
            cmd = "filter",
            args = '-tag "rq":">=' + min_read_quality + '"',
            runtime_attr_override = disable_preemption_runtime_attrs
    }

    call Utils.Bamtools as t_118_GetCcsReclaimedReadsWithCorrectedUmis {
        input:
            bamfile =t_116_CorrectUmisWithSetCover.corrected_umi_reads,
            prefix = SM + "_annotated_array_elements_for_quant_with_gene_names.corrected_umis.CCS_Reclaimed",
            cmd = "filter",
            args = '-tag "rq":"<' + min_read_quality + '"',
            runtime_attr_override = disable_preemption_runtime_attrs
    }

    call UMI_TOOLS.Run_Group as t_119_UMIToolsGroup {
        input:
            aligned_transcriptome_reads = t_116_CorrectUmisWithSetCover.corrected_umi_reads,
            aligned_transcriptome_reads_index = t_116_CorrectUmisWithSetCover.corrected_umi_reads_index,
            do_per_cell = true,
            prefix = SM + "_annotated_array_elements_with_gene_names_with_umi_tools_group_correction"
    }

    # Create CCS count matrix and anndata:
    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_120_CreateCCSCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_117_GetCcsCorrectedReadsWithCorrectedUmis.bam_out,
            tx_equivalence_class_assignments = t_114_CombineEqClassFiles.combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = SM + "_ccs_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_121_CreateCCSCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_120_CreateCCSCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = t_101_ST2_Quant.st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = t_114_CombineEqClassFiles.combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = t_114_CombineEqClassFiles.combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = t_114_CombineEqClassFiles.combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = t_114_CombineEqClassFiles.combined_gene_eq_class_defs,
            prefix = SM + "_ccs_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: 64}
    }

    # Create CLR count matrix and anndata:
    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_122_CreateCLRCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_118_GetCcsReclaimedReadsWithCorrectedUmis.bam_out,
            tx_equivalence_class_assignments = t_114_CombineEqClassFiles.combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = SM + "_clr_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_123_CreateCLRCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_122_CreateCLRCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = t_101_ST2_Quant.st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = t_114_CombineEqClassFiles.combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = t_114_CombineEqClassFiles.combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = t_114_CombineEqClassFiles.combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = t_114_CombineEqClassFiles.combined_gene_eq_class_defs,
            prefix = SM + "_clr_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: 64}
    }

    # Create overall count matrix and anndata:
    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_124_CreateOverallCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_116_CorrectUmisWithSetCover.corrected_umi_reads,
            tx_equivalence_class_assignments = t_114_CombineEqClassFiles.combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = SM + "_overall_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_125_CreateOverallCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_124_CreateOverallCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = t_101_ST2_Quant.st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = t_114_CombineEqClassFiles.combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = t_114_CombineEqClassFiles.combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = t_114_CombineEqClassFiles.combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = t_114_CombineEqClassFiles.combined_gene_eq_class_defs,
            prefix = SM + "_overall_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: 64}
    }


    #################################################
    #     ___   ____      __  __  __ _____ _____ ____  ___ ____ ____
    #    / _ \ / ___|    / / |  \/  | ____|_   _|  _ \|_ _/ ___/ ___|
    #   | | | | |       / /  | |\/| |  _|   | | | |_) || | |   \___ \
    #   | |_| | |___   / /   | |  | | |___  | | |  _ < | | |___ ___) |
    #    \__\_\\____| /_/    |_|  |_|_____| |_| |_| \_\___\____|____/
    #
    #################################################

    call AM.SamtoolsStats as t_126_BaselineArrayElementStats {
        input:
            bam = t_81_MergeAllArrayElementsNonTruncated.merged_bam
    }

    call AM.SamtoolsStats as t_127_AlignedArrayElementStats {
        input:
            bam = t_95_MergeAllAlignedArrayElementsNonTruncated.merged_bam
    }

    call AM.SamtoolsStats as t_128_AlignedFilteredArrayElementStats {
        input:
            bam = t_100_MergeAllAlignedAndFilteredArrayElements.merged_bam
    }

    call AM.SamtoolsStats as t_129_AlignedAnnotatedArrayElementsForQuantStats {
        input:
            bam = t_116_CorrectUmisWithSetCover.corrected_umi_reads
    }

    call LONGBOW.AggregateCorrectLogStats as t_130_AggregateLongbowCorrectStats {
        input:
            longbow_correct_log_files = flatten([flatten(t_45_LongbowCorrectCcsReclaimedArrayElementCBCs.log), flatten(t_28_LongbowCorrectCCSCorrectedArrayElementCBCs.log)]),
            out_name = SM + "_longbow_correct_stats.txt"
    }

    # Get stats on CCS reads:
    call LONGBOW.Stats as t_131_CCS_longbow_stats {
        input:
            reads = t_59_MergeCCSLongbowAnnotatedArrayReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_CCS_Corrected",
    }

    # Get stats on Reclaimable reads:
    call LONGBOW.Stats as t_132_Reclaimable_longbow_stats {
        input:
            reads = t_61_MergeCCSReclaimableLongbowAnnotatedArrayReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_CCS_Reclaimable",
    }

    # Get stats on Reclaimed reads:
    call LONGBOW.Stats as t_133_Reclaimed_longbow_stats {
        input:
            reads = t_67_MergeCCSReclaimedArrayReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_CCS_Reclaimed",
    }

    # Get stats on All Passing reads (overall stats):
    call LONGBOW.Stats as t_134_Passed_longbow_stats {
        input:
            reads = t_71_MergeLongbowPassedReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_All_Longbow_Passed",
    }

    # Get stats on All Failed reads (overall stats):
    call LONGBOW.Stats as t_135_Failed_longbow_stats {
        input:
            reads = t_73_MergeLongbowFailedReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_All_Longbow_Failed",
    }

    # Get stats on All reads (overall stats):
    call LONGBOW.Stats as t_136_Overall_longbow_stats {
        input:
            reads = t_75_MergeAllLongbowAnnotatedReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_Overall",
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

#    File keyfile = t_125_CreateOverallCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad

    # This seems to take longer to get to:
    File keyfile = t_119_UMIToolsGroup.output_tsv

    String base_out_dir = outdir + "/" + DIR + out_dir_suffix + "/" + t_01_WdlExecutionStartTimestamp.timestamp_string
    String stats_out_dir = base_out_dir + "/stats"
    String array_element_dir = base_out_dir + "/annotated_array_elements"
    String intermediate_reads_dir = base_out_dir + "/intermediate_reads"

    String meta_files_dir = base_out_dir + "/meta_files"

    String intermediate_array_reads_dir = intermediate_reads_dir + "/array_reads"
    String intermediate_array_elements_dir = intermediate_reads_dir + "/array_elements"

    String quant_dir = base_out_dir + "/quant"

    ##############################################################################################################
    # Finalize gene / tx assignments:
    call FF.FinalizeToDir as t_137_FinalizeEqClasses {
        input:
            files = [
                t_114_CombineEqClassFiles.combined_gene_eq_class_defs,
                t_114_CombineEqClassFiles.combined_gene_eq_class_assignments,
                t_114_CombineEqClassFiles.combined_tx_eq_class_defs,
                t_114_CombineEqClassFiles.combined_tx_eq_class_assignments,
            ],
            outdir = quant_dir + "/eqivalence_classes",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_138_FinalizeUmiToolsOutputs {
        input:
            files = [
                t_119_UMIToolsGroup.output_bam,
                t_119_UMIToolsGroup.output_tsv,
            ],
            outdir = quant_dir + "/UMITools",
            keyfile = keyfile
    }

    # CCS:
    call FF.FinalizeToDir as t_139_FinalizeCCSTxAndGeneAssignments {
        input:
            files = [
                t_120_CreateCCSCountMatrixFromAnnotatedBam.count_matrix,
                t_121_CreateCCSCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad,
            ],
            outdir = quant_dir + "/CCS",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_140_FinalizeCCSRawQuantPickles {
        input:
            files = t_121_CreateCCSCountMatrixAnndataFromEqClasses.pickles,
            outdir = quant_dir + "/CCS",
            keyfile = keyfile
    }

    # CLR:
    call FF.FinalizeToDir as t_141_FinalizeCLRTxAndGeneAssignments {
        input:
            files = [
                t_122_CreateCLRCountMatrixFromAnnotatedBam.count_matrix,
                t_123_CreateCLRCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad,
            ],
            outdir = quant_dir + "/CLR",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_142_FinalizeCLRRawQuantPickles {
        input:
            files = t_123_CreateCLRCountMatrixAnndataFromEqClasses.pickles,
            outdir = quant_dir + "/CLR",
            keyfile = keyfile
    }

    # Overall:
    call FF.FinalizeToDir as t_143_FinalizeOverallTxAndGeneAssignments {
        input:
            files = [
                t_124_CreateOverallCountMatrixFromAnnotatedBam.count_matrix,
                t_125_CreateOverallCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad,
            ],
            outdir = quant_dir + "/Overall",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_144_FinalizeOverallRawQuantPickles {
        input:
            files = t_125_CreateOverallCountMatrixAnndataFromEqClasses.pickles,
            outdir = quant_dir + "/Overall",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_145_FinalizeRefAndSt2Comparisons {
        input:
            files = [
                t_107_GffCompareStringtie2toGencode.refmap,
                t_107_GffCompareStringtie2toGencode.tmap,
                t_107_GffCompareStringtie2toGencode.tracking,
                t_107_GffCompareStringtie2toGencode.loci,
                t_107_GffCompareStringtie2toGencode.annotated_gtf,
                t_107_GffCompareStringtie2toGencode.stats,
                t_107_GffCompareStringtie2toGencode.log,

                t_108_GffCompareGencodetoStringtie2.refmap,
                t_108_GffCompareGencodetoStringtie2.tmap,
                t_108_GffCompareGencodetoStringtie2.tracking,
                t_108_GffCompareGencodetoStringtie2.loci,
                t_108_GffCompareGencodetoStringtie2.annotated_gtf,
                t_108_GffCompareGencodetoStringtie2.stats,
                t_108_GffCompareGencodetoStringtie2.log,
            ],
            outdir = quant_dir + "/gencode_and_stringtie2",
            keyfile = keyfile
    }

    # Finalize gene / tx assignment by contig:
    # NOTE: According to the scatter/gather documentation in the WDL spec, this will work correctly
    #       (https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#scatter--gather)
    scatter (i in range(length(t_109_SplitArrayElementsByContig.contig_bams))) {
        String contig = t_109_SplitArrayElementsByContig.contig_names[i]

        call FF.FinalizeToDir as t_146_FinalizeTxAndGeneAssignmentsByContig {
            input:
                files = [
                    t_111_GffCompareStringtie2toMasSeqReads.refmap[i],
                    t_111_GffCompareStringtie2toMasSeqReads.tmap[i],
                    t_111_GffCompareStringtie2toMasSeqReads.tracking[i],
                    t_111_GffCompareStringtie2toMasSeqReads.loci[i],
                    t_111_GffCompareStringtie2toMasSeqReads.annotated_gtf[i],
                    t_111_GffCompareStringtie2toMasSeqReads.stats[i],
                    t_111_GffCompareStringtie2toMasSeqReads.log[i],

                    t_112_GffCompareGencodetoMasSeqReads.refmap[i],
                    t_112_GffCompareGencodetoMasSeqReads.tmap[i],
                    t_112_GffCompareGencodetoMasSeqReads.tracking[i],
                    t_112_GffCompareGencodetoMasSeqReads.loci[i],
                    t_112_GffCompareGencodetoMasSeqReads.annotated_gtf[i],
                    t_112_GffCompareGencodetoMasSeqReads.stats[i],
                    t_112_GffCompareGencodetoMasSeqReads.log[i],

                    t_113_QuantifyGffComparison.gene_assignments_file[i],
                    t_113_QuantifyGffComparison.gene_eq_class_labels_file[i],
                    t_113_QuantifyGffComparison.tx_equivalence_class_labels_file[i],
                    t_113_QuantifyGffComparison.tx_equivalence_class_file[i],
                    t_113_QuantifyGffComparison.graph_gpickle[i],
                ],
                outdir = quant_dir + "/by_contig/" + contig,
                keyfile = keyfile
        }
    }

    ##############################################################################################################
    # Finalize annotated, aligned array elements:
    call FF.FinalizeToDir as t_147_FinalizeIntermediateAnnotatedArrayElements {
        input:
            files = [
                t_77_MergeCCSArrayElements.merged_bam,
                t_77_MergeCCSArrayElements.merged_bai,
                t_79_MergeCCSArrayElementsNonTruncated.merged_bam,
                t_79_MergeCCSArrayElementsNonTruncated.merged_bai,
                t_82_MergeCCSArrayElementsUmiPadded.merged_bam,
                t_82_MergeCCSArrayElementsUmiPadded.merged_bai,
                t_83_MergeCCSArrayElementsUmiCbcPadded.merged_bam,
                t_83_MergeCCSArrayElementsUmiCbcPadded.merged_bai,
                t_84_MergeCCSArrayElementsUmiCbcPaddedCbcCorrected.merged_bam,
                t_84_MergeCCSArrayElementsUmiCbcPaddedCbcCorrected.merged_bai,
                t_85_MergeLongbowPaddedCBCCorrectedCCSArrayElements.merged_bam,
                t_85_MergeLongbowPaddedCBCCorrectedCCSArrayElements.merged_bai,
                t_86_MergeLongbowPaddedCBCUncorrectableCCSArrayElements.merged_bam,
                t_86_MergeLongbowPaddedCBCUncorrectableCCSArrayElements.merged_bai,

                t_78_MergeCCSReclaimedArrayElements.merged_bam,
                t_78_MergeCCSReclaimedArrayElements.merged_bai,
                t_80_MergeCCSReclaimedArrayElementsNonTruncated.merged_bam,
                t_80_MergeCCSReclaimedArrayElementsNonTruncated.merged_bai,
                t_88_MergeCCSReclaimedArrayElementsUmiPadded.merged_bam,
                t_88_MergeCCSReclaimedArrayElementsUmiPadded.merged_bai,
                t_89_MergeCCSReclaimedArrayElementsUmiCbcPadded.merged_bam,
                t_89_MergeCCSReclaimedArrayElementsUmiCbcPadded.merged_bai,
                t_90_MergeCCSReclaimedArrayElementsUmiCbcPaddedCbcCorrected.merged_bam,
                t_90_MergeCCSReclaimedArrayElementsUmiCbcPaddedCbcCorrected.merged_bai,
                t_91_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElements.merged_bam,
                t_91_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElements.merged_bai,

                t_81_MergeAllArrayElementsNonTruncated.merged_bam,
                t_81_MergeAllArrayElementsNonTruncated.merged_bai,

                t_95_MergeAllAlignedArrayElementsNonTruncated.merged_bam,
                t_95_MergeAllAlignedArrayElementsNonTruncated.merged_bai,

                t_100_MergeAllAlignedAndFilteredArrayElements.merged_bam,
                t_100_MergeAllAlignedAndFilteredArrayElements.merged_bai,

                t_85_MergeLongbowPaddedCBCCorrectedCCSArrayElements.merged_bam,
                t_85_MergeLongbowPaddedCBCCorrectedCCSArrayElements.merged_bai,
                t_86_MergeLongbowPaddedCBCUncorrectableCCSArrayElements.merged_bam,
                t_86_MergeLongbowPaddedCBCUncorrectableCCSArrayElements.merged_bai,
                t_87_MergeLongbowPaddedCBCCorrectedExtractedCCSArrayElements.merged_bam,
                t_87_MergeLongbowPaddedCBCCorrectedExtractedCCSArrayElements.merged_bai,
                t_90_MergeCCSReclaimedArrayElementsUmiCbcPaddedCbcCorrected.merged_bam,
                t_90_MergeCCSReclaimedArrayElementsUmiCbcPaddedCbcCorrected.merged_bai,
                t_91_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElements.merged_bam,
                t_91_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElements.merged_bai,
                t_92_MergeLongbowPaddedCBCCorrectedExtractedCCSReclaimedArrayElements.merged_bam,
                t_92_MergeLongbowPaddedCBCCorrectedExtractedCCSReclaimedArrayElements.merged_bai,

                t_104_RestoreCcsOriginalReadNames.bam_out,
                t_105_RestoreClrOriginalReadNames.bam_out,

                t_106_MergeAllAnnotatedArrayElementsWithOriginalNames.merged_bam,
                t_106_MergeAllAnnotatedArrayElementsWithOriginalNames.merged_bai,

                t_116_CorrectUmisWithSetCover.uncorrected_umi_reads
            ],
            outdir = intermediate_array_elements_dir,
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_148_FinalizeAnnotatedArrayElements {
        input:
            files = [
                t_116_CorrectUmisWithSetCover.corrected_umi_reads,
                t_116_CorrectUmisWithSetCover.corrected_umi_reads_index,

                t_117_GetCcsCorrectedReadsWithCorrectedUmis.bam_out,
                t_118_GetCcsReclaimedReadsWithCorrectedUmis.bam_out,

                t_100_MergeAllAlignedAndFilteredArrayElements.merged_bam,
                t_100_MergeAllAlignedAndFilteredArrayElements.merged_bai
            ],
            outdir = array_element_dir,
            keyfile = keyfile
    }

    ##############################################################################################################
    # Finalize meta files:
    call FF.FinalizeToDir as t_149_FinalizeMeta {
        input:
            files = [
                cell_barcode_whitelist,
                t_96_MergeAllCCSBarcodeConfShards.merged_file,
                t_97_MergeAllCCSReclaimedBarcodeConfShards.merged_file
            ],
            outdir = meta_files_dir,
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_150_FinalizeCCSCBCcorrectionLogsToMeta {
        input:
            files = flatten(t_28_LongbowCorrectCCSCorrectedArrayElementCBCs.log),
            outdir = meta_files_dir + "/" + "ccs_cbc_correction_logs",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_151_FinalizeCCSRejectedCBCcorrectionLogsToMeta {
        input:
            files = flatten(t_45_LongbowCorrectCcsReclaimedArrayElementCBCs.log),
            outdir = meta_files_dir + "/" + "ccs_rejected_cbc_correction_logs",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_152_FinalizeCCSUmiAdjustmentLogs {
        input:
            files = flatten(t_29_AdjustCCSUMIs.log),
            outdir = meta_files_dir + "/" + "umi_adjustment_logs_ccs",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_153_FinalizeCCSReclaimedUmiAdjustmentLogs {
        input:
            files = flatten(t_46_AdjustCcsReclaimedUMIs.log),
            outdir = meta_files_dir + "/" + "umi_adjustment_logs_ccs_reclaimed",
            keyfile = keyfile
    }

    ##############################################################################################################
    # Finalize the discovered transcriptome:
    if ( !is_SIRV_data ) {
        call FF.FinalizeToDir as t_154_FinalizeDiscoveredTranscriptome {
            input:
                files = [
                    t_101_ST2_Quant.st_gtf,
                    t_102_ST2_ExtractTranscriptSequences.transcripts_fa,
                    t_102_ST2_ExtractTranscriptSequences.transcripts_fai,
                    t_102_ST2_ExtractTranscriptSequences.transcripts_dict,
                    t_103_ST2_CompareTranscriptomes.annotated_gtf,
                    t_103_ST2_CompareTranscriptomes.loci,
                    t_103_ST2_CompareTranscriptomes.stats,
                    t_103_ST2_CompareTranscriptomes.tracking,
                    t_103_ST2_CompareTranscriptomes.refmap,
                    t_103_ST2_CompareTranscriptomes.tmap,
                ],
                outdir = base_out_dir + "/discovered_transcriptome",
                keyfile = keyfile
        }
    }
    ##############################################################################################################
    # Finalize the intermediate reads files (from raw CCS corrected reads through split array elements)
    call FF.FinalizeToDir as t_155_FinalizeArrayReads {
        input:
            files = [

                t_59_MergeCCSLongbowAnnotatedArrayReads.merged_bam,
                t_59_MergeCCSLongbowAnnotatedArrayReads.merged_bai,
                t_60_PbIndexMergedCCSLongbowAnnotatedArrayReads.pbindex,

                t_61_MergeCCSReclaimableLongbowAnnotatedArrayReads.merged_bam,
                t_61_MergeCCSReclaimableLongbowAnnotatedArrayReads.merged_bai,
                t_62_PbIndexMergedCCSReclaimableLongbowAnnotatedArrayReads.pbindex,

                t_71_MergeLongbowPassedReads.merged_bam,
                t_71_MergeLongbowPassedReads.merged_bai,
                t_72_PbIndexMergedLongbowPassingReads.pbindex,

                t_73_MergeLongbowFailedReads.merged_bam,
                t_73_MergeLongbowFailedReads.merged_bai,
                t_74_PbIndexMergedLongbowFailedReads.pbindex,

                t_75_MergeAllLongbowAnnotatedReads.merged_bam,
                t_75_MergeAllLongbowAnnotatedReads.merged_bai,
                t_76_PbIndexMergedAllLongbowAnnotatedReads.pbindex,

                t_63_MergeCCSLongbowPassedArrayReads.merged_bam,
                t_63_MergeCCSLongbowPassedArrayReads.merged_bai,
                t_64_PbIndexMergedCCSLongbowPassedArrayReads.pbindex,

                t_65_MergeCCSLongbowFailedArrayReads.merged_bam,
                t_65_MergeCCSLongbowFailedArrayReads.merged_bai,
                t_66_PbIndexMergedCCSLongbowFailedArrayReads.pbindex,

                t_67_MergeCCSReclaimedArrayReads.merged_bam,
                t_67_MergeCCSReclaimedArrayReads.merged_bai,
                t_68_PbIndexMergedCCSReclaimedArrayReads.pbindex,

                t_69_MergeCCSUnreclaimableArrayReads.merged_bam,
                t_69_MergeCCSUnreclaimableArrayReads.merged_bai,
                t_70_PbIndexMergedCCSUnreclaimableReclaimedArrayReads.pbindex,

            ],
            outdir = intermediate_array_reads_dir,
            keyfile = keyfile
    }

    ##############################################################################################################
    # Finalize Stats:

    # Write out completion file so in the future we can be 100% sure that this run was good:
    call FF.FinalizeToDir as t_156_FinalizeHighLevelStats {
        input:
            files = [ t_11_FindCCSReport.ccs_report[0], t_130_AggregateLongbowCorrectStats.stats ],
            outdir = stats_out_dir,
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_157_FinalizeQuantArrayElementStats {
        input:
            files = [
                t_129_AlignedAnnotatedArrayElementsForQuantStats.raw_stats,
                t_129_AlignedAnnotatedArrayElementsForQuantStats.summary_stats,
                t_129_AlignedAnnotatedArrayElementsForQuantStats.first_frag_qual,
                t_129_AlignedAnnotatedArrayElementsForQuantStats.last_frag_qual,
                t_129_AlignedAnnotatedArrayElementsForQuantStats.first_frag_gc_content,
                t_129_AlignedAnnotatedArrayElementsForQuantStats.last_frag_gc_content,
                t_129_AlignedAnnotatedArrayElementsForQuantStats.acgt_content_per_cycle,
                t_129_AlignedAnnotatedArrayElementsForQuantStats.insert_size,
                t_129_AlignedAnnotatedArrayElementsForQuantStats.read_length_dist,
                t_129_AlignedAnnotatedArrayElementsForQuantStats.indel_distribution,
                t_129_AlignedAnnotatedArrayElementsForQuantStats.indels_per_cycle,
                t_129_AlignedAnnotatedArrayElementsForQuantStats.coverage_distribution,
                t_129_AlignedAnnotatedArrayElementsForQuantStats.gc_depth,
            ],
            outdir = stats_out_dir + "/array_elements_for_quant/",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_158_FinalizeTxomeDiscoveryArrayElementStats {
        input:
            files = [
                t_128_AlignedFilteredArrayElementStats.raw_stats,
                t_128_AlignedFilteredArrayElementStats.summary_stats,
                t_128_AlignedFilteredArrayElementStats.first_frag_qual,
                t_128_AlignedFilteredArrayElementStats.last_frag_qual,
                t_128_AlignedFilteredArrayElementStats.first_frag_gc_content,
                t_128_AlignedFilteredArrayElementStats.last_frag_gc_content,
                t_128_AlignedFilteredArrayElementStats.acgt_content_per_cycle,
                t_128_AlignedFilteredArrayElementStats.insert_size,
                t_128_AlignedFilteredArrayElementStats.read_length_dist,
                t_128_AlignedFilteredArrayElementStats.indel_distribution,
                t_128_AlignedFilteredArrayElementStats.indels_per_cycle,
                t_128_AlignedFilteredArrayElementStats.coverage_distribution,
                t_128_AlignedFilteredArrayElementStats.gc_depth,
            ],
            outdir = stats_out_dir + "/array_elements_for_transcriptome_discovery/",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_159_FinalizeAlignedArrayElementStats {
        input:
            files = [
                t_127_AlignedArrayElementStats.raw_stats,
                t_127_AlignedArrayElementStats.summary_stats,
                t_127_AlignedArrayElementStats.first_frag_qual,
                t_127_AlignedArrayElementStats.last_frag_qual,
                t_127_AlignedArrayElementStats.first_frag_gc_content,
                t_127_AlignedArrayElementStats.last_frag_gc_content,
                t_127_AlignedArrayElementStats.acgt_content_per_cycle,
                t_127_AlignedArrayElementStats.insert_size,
                t_127_AlignedArrayElementStats.read_length_dist,
                t_127_AlignedArrayElementStats.indel_distribution,
                t_127_AlignedArrayElementStats.indels_per_cycle,
                t_127_AlignedArrayElementStats.coverage_distribution,
                t_127_AlignedArrayElementStats.gc_depth,
            ],
            outdir = stats_out_dir + "/aligned_array_elements/",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_160_FinalizeBaselineArrayElementStats {
        input:
            files = [
                t_126_BaselineArrayElementStats.raw_stats,
                t_126_BaselineArrayElementStats.summary_stats,
                t_126_BaselineArrayElementStats.first_frag_qual,
                t_126_BaselineArrayElementStats.last_frag_qual,
                t_126_BaselineArrayElementStats.first_frag_gc_content,
                t_126_BaselineArrayElementStats.last_frag_gc_content,
                t_126_BaselineArrayElementStats.acgt_content_per_cycle,
                t_126_BaselineArrayElementStats.insert_size,
                t_126_BaselineArrayElementStats.read_length_dist,
                t_126_BaselineArrayElementStats.indel_distribution,
                t_126_BaselineArrayElementStats.indels_per_cycle,
                t_126_BaselineArrayElementStats.coverage_distribution,
                t_126_BaselineArrayElementStats.gc_depth,
            ],
            outdir = stats_out_dir + "/baseline_array_elements/",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_161_FinalizeCCSLongbowStats {
        input:
            files = [
                t_131_CCS_longbow_stats.summary_stats,
                t_131_CCS_longbow_stats.array_length_counts_plot_png,
                t_131_CCS_longbow_stats.array_length_counts_plot_svg,
                t_131_CCS_longbow_stats.ligation_heatmap_nn_png,
                t_131_CCS_longbow_stats.ligation_heatmap_nn_svg,
                t_131_CCS_longbow_stats.ligation_heatmap_png,
                t_131_CCS_longbow_stats.ligation_heatmap_svg,
                t_131_CCS_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_131_CCS_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_131_CCS_longbow_stats.ligation_heatmap_reduced_png,
                t_131_CCS_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/CCS_Corrected/",
            keyfile = keyfile
    }
    call FF.FinalizeToDir as t_162_FinalizeReclaimableLongbowStats {
        input:
            files = [
                t_132_Reclaimable_longbow_stats.summary_stats,
                t_132_Reclaimable_longbow_stats.array_length_counts_plot_png,
                t_132_Reclaimable_longbow_stats.array_length_counts_plot_svg,
                t_132_Reclaimable_longbow_stats.ligation_heatmap_nn_png,
                t_132_Reclaimable_longbow_stats.ligation_heatmap_nn_svg,
                t_132_Reclaimable_longbow_stats.ligation_heatmap_png,
                t_132_Reclaimable_longbow_stats.ligation_heatmap_svg,
                t_132_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_132_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_132_Reclaimable_longbow_stats.ligation_heatmap_reduced_png,
                t_132_Reclaimable_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/CCS_Reclaimable/",
            keyfile = keyfile
    }
    call FF.FinalizeToDir as t_163_FinalizeReclaimedLongbowStats {
        input:
            files = [
                t_133_Reclaimed_longbow_stats.summary_stats,
                t_133_Reclaimed_longbow_stats.array_length_counts_plot_png,
                t_133_Reclaimed_longbow_stats.array_length_counts_plot_svg,
                t_133_Reclaimed_longbow_stats.ligation_heatmap_nn_png,
                t_133_Reclaimed_longbow_stats.ligation_heatmap_nn_svg,
                t_133_Reclaimed_longbow_stats.ligation_heatmap_png,
                t_133_Reclaimed_longbow_stats.ligation_heatmap_svg,
                t_133_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_133_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_133_Reclaimed_longbow_stats.ligation_heatmap_reduced_png,
                t_133_Reclaimed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/CCS_Reclaimed/",
            keyfile = keyfile
    }
    call FF.FinalizeToDir as t_164_FinalizeOverallLongbowStats {
        input:
            files = [
                t_136_Overall_longbow_stats.summary_stats,
                t_136_Overall_longbow_stats.array_length_counts_plot_png,
                t_136_Overall_longbow_stats.array_length_counts_plot_svg,
                t_136_Overall_longbow_stats.ligation_heatmap_nn_png,
                t_136_Overall_longbow_stats.ligation_heatmap_nn_svg,
                t_136_Overall_longbow_stats.ligation_heatmap_png,
                t_136_Overall_longbow_stats.ligation_heatmap_svg,
                t_136_Overall_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_136_Overall_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_136_Overall_longbow_stats.ligation_heatmap_reduced_png,
                t_136_Overall_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/Overall/",
            keyfile = keyfile
    }
    call FF.FinalizeToDir as t_165_FinalizeAllPassedLongbowStats {
        input:
            files = [
                t_134_Passed_longbow_stats.summary_stats,
                t_134_Passed_longbow_stats.array_length_counts_plot_png,
                t_134_Passed_longbow_stats.array_length_counts_plot_svg,
                t_134_Passed_longbow_stats.ligation_heatmap_nn_png,
                t_134_Passed_longbow_stats.ligation_heatmap_nn_svg,
                t_134_Passed_longbow_stats.ligation_heatmap_png,
                t_134_Passed_longbow_stats.ligation_heatmap_svg,
                t_134_Passed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_134_Passed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_134_Passed_longbow_stats.ligation_heatmap_reduced_png,
                t_134_Passed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/All_Longbow_Passed/",
            keyfile = keyfile
    }
    call FF.FinalizeToDir as t_166_FinalizeAllPassedLongbowStats {
        input:
            files = [
                t_135_Failed_longbow_stats.summary_stats,
                t_135_Failed_longbow_stats.array_length_counts_plot_png,
                t_135_Failed_longbow_stats.array_length_counts_plot_svg,
                t_135_Failed_longbow_stats.ligation_heatmap_nn_png,
                t_135_Failed_longbow_stats.ligation_heatmap_nn_svg,
                t_135_Failed_longbow_stats.ligation_heatmap_png,
                t_135_Failed_longbow_stats.ligation_heatmap_svg,
                t_135_Failed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_135_Failed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_135_Failed_longbow_stats.ligation_heatmap_reduced_png,
                t_135_Failed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/All_Longbow_Failed/",
            keyfile = keyfile
    }

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good:
    call FF.WriteCompletionFile as t_167_WriteCompletionFile {
        input:
            outdir = base_out_dir + "/",
            keyfile = keyfile
    }
}
