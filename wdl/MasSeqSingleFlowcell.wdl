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

workflow MasSeqSingleFlowcell {

    meta {
        description : "This workflow is designed to process data from the MASSeq v2 protocol and produce aligned reads that are ready for downstream analysis (e.g. transcript isoform identification).  It takes in a raw PacBio run folder location on GCS and produces a folder containing the aligned reads and other processed data."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String gcs_input_dir
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/MasSeqSingleFlowcell"

        File head_adapter_fasta = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/10x_adapter.fasta"
        File tail_adapter_fasta = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/tso_adapter.fasta"
        File ten_x_cell_barcode_whitelist = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/737K-august-2016.txt"

        # NOTE: Reference for un-split CCS reads:
        File ref_fasta =  "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
        File ref_fasta_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
        File ref_fasta_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"

        # NOTE: Reference for array elements:
        File transcriptome_ref_fasta =  "gs://broad-dsde-methods-long-reads-public/resources/gencode_v37/gencode.v37.pc_transcripts.fa"
        File transcriptome_ref_fasta_index = "gs://broad-dsde-methods-long-reads-public/resources/gencode_v37/gencode.v37.pc_transcripts.fa.fai"
        File transcriptome_ref_fasta_dict = "gs://broad-dsde-methods-long-reads-public/resources/gencode_v37/gencode.v37.pc_transcripts.dict"

        File genome_annotation_gtf = "gs://broad-dsde-methods-long-reads-public/resources/gencode_v37/gencode.v37.primary_assembly.annotation.gtf"

        File jupyter_template_static = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/MAS-seq_QC_report_template-static.ipynb"
        File workflow_dot_file = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/MasSeqSingleFlowcell.dot"

        File intervals_of_interest = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/gencode.37.TCR_intervals.tsv"
        String interval_overlap_name = "is_tcr_overlapping"

        File? illumina_barcoded_bam

        # Default here is 0 because ccs uncorrected reads all seem to have RQ = -1.
        # All pathologically long reads also have RQ = -1.
        # This way we preserve the vast majority of the data, even if it has low quality.
        # We can filter it out at later steps.
        Float min_read_quality = 0.0
        Int max_reclamation_length = 60000

        Boolean is_SIRV_data = false
        String mas_seq_model = "mas15"

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

    call PB.GetRunInfo as t_06_GetRunInfo { input: subread_bam = reads_bam }

    String SM  = select_first([sample_name, t_06_GetRunInfo.run_info["SM"]])
    String PL  = "PACBIO"
    String PU  = t_06_GetRunInfo.run_info["PU"]
    String DT  = t_06_GetRunInfo.run_info["DT"]
    String ID  = PU
    String DS  = t_06_GetRunInfo.run_info["DS"]
    String DIR = SM + "." + ID

    String RG_subreads  = "@RG\\tID:~{ID}.subreads\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
    String RG_consensus = "@RG\\tID:~{ID}.consensus\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
    String RG_array_elements = "@RG\\tID:~{ID}.array_elements\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

    # Check to see if we need to annotate our reads:
    call LONGBOW.CheckForAnnotatedArrayReads as t_07_CheckForAnnotatedReads {
        input:
            bam = reads_bam
    }

    File read_pbi = sub(reads_bam, ".bam$", ".bam.pbi")
    call PB.ShardLongReads as t_08_ShardLongReads {
        input:
            unaligned_bam = reads_bam,
            unaligned_pbi = read_pbi,
            prefix = SM + "_shard",
            num_shards = 300,
    }

    scatter (sharded_reads in t_08_ShardLongReads.unmapped_shards) {

        ## No more preemption on this sharding - takes too long otherwise.
        RuntimeAttr disable_preemption_runtime_attrs = object {
            preemptible_tries: 0
        }

        String fbmrq_prefix = basename(sharded_reads, ".bam")

        # Filter out the kinetics tags from PB files:
        call PB.RemoveKineticsTags as t_09_RemoveKineticsTags {
            input:
                bam = sharded_reads,
                prefix = SM + "_kinetics_removed"
        }

        # Handle setting up the things that we need for further processing of CCS-only reads:
        call PB.FindCCSReport as t_10_FindCCSReport {
            input:
                gcs_input_dir = gcs_input_dir
        }

        # 1 - filter the reads by the minimum read quality:
        call Utils.Bamtools as t_11_FilterS2EByMinReadQuality {
            input:
                bamfile = t_09_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_good_reads",
                cmd = "filter",
                args = '-tag "rq":">=' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # 1.5 - Get the "rejected" reads:
        call Utils.Bamtools as t_12_GetS2ERCcsRejectedReads {
            input:
                bamfile = t_09_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_rejected_reads",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # 2 - Get reads we can reclaim:
        call Utils.Bamtools as t_13_ExtractS2ECcsReclaimableReads {
            input:
                bamfile = t_09_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_reads_for_ccs_reclamation",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '" -length "<=' + max_reclamation_length + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        if ( ! t_07_CheckForAnnotatedReads.bam_has_annotations ) {
            # 3: Longbow annotate ccs reads
            call LONGBOW.Annotate as t_14_AnnotateS2ECCSReads {
                input:
                    reads = t_11_FilterS2EByMinReadQuality.bam_out,
                    model = mas_seq_model
            }
            # 4: Longbow annotate reclaimable reads
            call LONGBOW.Annotate as t_15_AnnotateS2EReclaimableReads {
                input:
                    reads = t_13_ExtractS2ECcsReclaimableReads.bam_out,
                    model = mas_seq_model
            }
        }

        File annotated_S2E_ccs_file = if t_07_CheckForAnnotatedReads.bam_has_annotations then t_11_FilterS2EByMinReadQuality.bam_out else select_first([t_14_AnnotateS2ECCSReads.annotated_bam])
        File annotated_S2E_reclaimable_file = if t_07_CheckForAnnotatedReads.bam_has_annotations then t_13_ExtractS2ECcsReclaimableReads.bam_out else select_first([t_15_AnnotateS2EReclaimableReads.annotated_bam])

        # 5: Longbow filter ccs annotated reads
        call LONGBOW.Filter as t_16_FilterS2ECCSReads {
            input:
                bam = annotated_S2E_ccs_file,
                prefix = SM + "_subshard",
                model = mas_seq_model
        }

        # 6: Longbow filter ccs reclaimable reads
        call LONGBOW.Filter as t_17_FilterS2EReclaimableReads {
            input:
                bam = annotated_S2E_reclaimable_file,
                prefix = SM + "_subshard",
                model = mas_seq_model
        }

        # 7: Merge reclaimed and ccs longbow filtered reads
        call Utils.MergeBams as t_18_MergeLongbowS2EPassedReads {
            input:
                bams = [t_16_FilterS2ECCSReads.passed_reads, t_17_FilterS2EReclaimableReads.passed_reads],
                prefix = SM + "_LongbowFilter_Passed_1"
        }
        call Utils.MergeBams as t_19_MergeLongbowS2EFailedReads {
            input:
                bams = [t_16_FilterS2ECCSReads.failed_reads, t_17_FilterS2EReclaimableReads.failed_reads],
                prefix = SM + "_LongbowFilter_Failed_1"
        }

        # 8: PBIndex reads
        call PB.PBIndex as t_20_PbIndexS2ELongbowPassedReads {
            input:
                bam = t_18_MergeLongbowS2EPassedReads.merged_bam
        }

        # 9: Get CCS Reclaimed array elements for further study:
        call PB.PBIndex as t_21_PbIndexS2ECcsReclaimedReads {
            input:
                bam = t_17_FilterS2EReclaimableReads.passed_reads
        }
        call PB.ShardLongReads as t_22_ShardS2ECcsReclaimedReads {
            input:
                unaligned_bam = t_17_FilterS2EReclaimableReads.passed_reads,
                unaligned_pbi = t_21_PbIndexS2ECcsReclaimedReads.pbindex,
                prefix = SM + "_ccs_reclaimed_reads_subshard",
                num_shards = 10,
        }
        scatter (s2e_ccs_reclaimed_shard in t_22_ShardS2ECcsReclaimedReads.unmapped_shards) {
            call LONGBOW.Segment as t_23_SegmentS2ECcsReclaimedReads {
                input:
                    annotated_reads = s2e_ccs_reclaimed_shard,
                    prefix = SM + "_ccs_reclaimed_array_elements_subshard",
                    model = mas_seq_model
            }
        }
        call Utils.MergeBams as t_24_MergeS2ECcsReclaimedArrayElementSubshards {
            input:
                bams = t_23_SegmentS2ECcsReclaimedReads.segmented_bam,
                prefix = SM + "_ccs_reclaimed_array_elements_shard"
        }

        # Shard these reads even wider so we can make sure we don't run out of memory:
        call PB.ShardLongReads as t_25_ShardCorrectedReads {
            input:
                unaligned_bam = t_18_MergeLongbowS2EPassedReads.merged_bam,
                unaligned_pbi = t_20_PbIndexS2ELongbowPassedReads.pbindex,
                prefix = SM + "_longbow_all_passed_subshard",
                num_shards = 10,
        }

        # Segment our arrays into individual array elements:
        scatter (corrected_shard in t_25_ShardCorrectedReads.unmapped_shards) {
            call LONGBOW.Segment as t_26_SegmentAnnotatedReads {
                input:
                    annotated_reads = corrected_shard,
                    model = mas_seq_model
            }
        }

        # Merge all outputs of Longbow Annotate / Segment:
        call Utils.MergeBams as t_27_MergeArrayElements_1 {
            input:
                bams = t_26_SegmentAnnotatedReads.segmented_bam,
                prefix = SM + "_ArrayElements_intermediate_1"
        }

        # The SIRV library prep is slightly different from the standard prep, so we have to account for it here:
        if (is_SIRV_data) {
            call TENX.TagSirvUmiPositionsFromLongbowAnnotatedArrayElement as t_28_TagSirvUmiPositionsFromLongbowAnnotatedArrayElements {
                input:
                    bam_file = t_27_MergeArrayElements_1.merged_bam,
                    prefix = SM + "_ArrayElements_SIRV_UMI_Extracted"
            }
        }
        if ( !is_SIRV_data ) {
            call TENX.AnnotateBarcodesAndUMIs as t_29_TenxAnnotateArrayElementsRaw {
                input:
                    bam_file = t_27_MergeArrayElements_1.merged_bam,
                    bam_index = t_27_MergeArrayElements_1.merged_bai,
                    head_adapter_fasta = head_adapter_fasta,
                    tail_adapter_fasta = tail_adapter_fasta,
                    whitelist_10x = ten_x_cell_barcode_whitelist,
                    read_end_length = 200,
                    poly_t_length = 31,
                    barcode_length = 16,
                    umi_length = 10,
                    raw_extract_only = true,
                    runtime_attr_override = fast_network_attrs
            }
        }
    }

    ############################################################################################################
    ############################################################################################################

    # Concatenate the TSV files with the barcode scores that we just created:
    if ( is_SIRV_data ) {
        call Utils.MergeBams as t_30_MergeSirvAnnotatedArrayElements {
            input:
                bams = select_all(t_28_TagSirvUmiPositionsFromLongbowAnnotatedArrayElements.output_bam),
                prefix = SM + "_SIRV_annotated_array_elements"
        }
    }
    if ( !is_SIRV_data ) {
        call Utils.MergeFiles as t_31_MergeUmiConfScoreTsvsForStarcode {
            input:
                files_to_merge = select_all(t_29_TenxAnnotateArrayElementsRaw.raw_starcode_counts)
        }

        # If we have our ilmn barcode file, we need to process it here:
        if (defined(illumina_barcoded_bam)) {
            call TENX.ExtractIlmnBarcodeConfScores as t_32_ExtractIlmnBarcodeConfScores {
                input:
                    bam_file = select_first([illumina_barcoded_bam]),
                    prefix = SM,
                    runtime_attr_override = fast_network_attrs
            }

            # Concatenate the TSV files with the barcode scores that we just created:
            call Utils.MergeFiles as t_33_GetMasterUmiConfScoreTsvForStarcode {
                input:
                    files_to_merge = [t_31_MergeUmiConfScoreTsvsForStarcode.merged_file, t_32_ExtractIlmnBarcodeConfScores.conf_score_tsv]
            }
        }
        File starcode_seeds = if (defined(illumina_barcoded_bam)) then select_first([t_33_GetMasterUmiConfScoreTsvForStarcode.merged_file]) else t_31_MergeUmiConfScoreTsvsForStarcode.merged_file

        # We have to consolidate our seeds into unique entries for starcode not to crash and burn:
        call TX_POST.MergeBarcodeCounts as t_34_ConsolidateBarcodeCountsForStarcode {
            input:
                barcode_count_tsv = starcode_seeds,
                prefix = SM + "_barcode_counts_for_starcode"
        }

        # Now we can correct our barcodes:
        scatter (tenx_annotated_bam in select_all(t_29_TenxAnnotateArrayElementsRaw.output_bam)) {
            call TENX.CorrectBarcodesWithStarcodeSeedCounts as t_35_CorrectBarcodesWithStarcodeSeedCounts {
                input:
                    bam_file = tenx_annotated_bam,
                    starcode_seeds_tsv = t_34_ConsolidateBarcodeCountsForStarcode.merged_counts,
                    whitelist_10x = ten_x_cell_barcode_whitelist,
                    prefix = SM + "_array_elements"
            }
        }

        # Merge the barcode corrected files here:
        call Utils.MergeBams as t_36_MergeAnnotatedArrayElements {
            input:
                bams = t_35_CorrectBarcodesWithStarcodeSeedCounts.output_bam,
                prefix = SM + "_annotated_array_elements"
        }
    }

    # Create an alias here that we can refer to in later steps regardless as to whether we have SIRV data or not
    # This `select_first` business is some sillyness to fix the conditional calls automatically converting the
    # output to `File?` instead of `File`
    File annotated_array_elements = if is_SIRV_data then select_first([t_30_MergeSirvAnnotatedArrayElements.merged_bam]) else select_first([t_36_MergeAnnotatedArrayElements.merged_bam])
    call PB.PBIndex as t_37_PbIndexAnnotatedArrayElements {
        input:
            bam = annotated_array_elements,
            runtime_attr_override = fast_network_attrs
    }

    ############################################################################################################
    ############################################################################################################
    # Now we can re-shard our array elements for extraction, alignment and labeling:

    call PB.ShardLongReads as t_38_ShardArrayElements {
        input:
            unaligned_bam = annotated_array_elements,
            unaligned_pbi = t_37_PbIndexAnnotatedArrayElements.pbindex,
            prefix = SM + "_ArrayElements_shard",
            num_shards = 300,
    }

    scatter (sharded_array_elements in t_38_ShardArrayElements.unmapped_shards) {

        # Grab only the coding regions of the annotated reads for our alignments:
        Int extract_start_offset = if is_SIRV_data then 8 else 26
        call LONGBOW.Extract as t_39_ExtractCodingRegionsFromArrayElements {
            input:
                bam = sharded_array_elements,
                start_offset = extract_start_offset,
                prefix = SM + "_ArrayElements_Coding_Regions_Only"
        }

        call AR.Minimap2 as t_40_AlignArrayElementsToGenome {
            input:
                reads      = [ t_39_ExtractCodingRegionsFromArrayElements.extracted_reads ],
                ref_fasta  = ref_fasta,
                map_preset = "splice:hq"
        }

        # We need to restore the annotations we created with the 10x tool to the aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_41_RestoreAnnotationsToGenomeAlignedBam {
            input:
                annotated_bam_file = t_39_ExtractCodingRegionsFromArrayElements.extracted_reads,
                aligned_bam_file = t_40_AlignArrayElementsToGenome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
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
    call Utils.MergeBams as t_42_MergeCCSRqFilteredReads { input: bams = t_11_FilterS2EByMinReadQuality.bam_out, prefix = SM + "_ccs_reads" }
    call Utils.MergeBams as t_43_MergeCCSRqRejectedReads { input: bams = t_12_GetS2ERCcsRejectedReads.bam_out, prefix = SM + "_ccs_rejected_reads" }
    call Utils.MergeBams as t_44_MergeAnnotatedCCSReads_S2e { input: bams = annotated_S2E_ccs_file, prefix = SM + "_ccs_reads_annotated" }
    call Utils.MergeBams as t_45_MergeLongbowPassedCCSReads_S2e { input: bams = t_16_FilterS2ECCSReads.passed_reads, prefix = SM + "_ccs_reads_annotated_longbow_passed" }
    call Utils.MergeBams as t_46_MergeLongbowFailedCCSReads_S2e { input: bams = t_16_FilterS2ECCSReads.failed_reads, prefix = SM + "_ccs_reads_annotated_longbow_failed" }

    # CCS Failed / Reclaimable:
    call Utils.MergeBams as t_47_MergeCCSReclaimableReads_S2e { input: bams = t_13_ExtractS2ECcsReclaimableReads.bam_out, prefix = SM + "_ccs_rejected_reclaimable" }
    call Utils.MergeBams as t_48_MergeCCSReclaimableAnnotatedReads_S2e { input: bams = annotated_S2E_reclaimable_file, prefix = SM + "_ccs_rejected_reclaimable_annotated" }
    call Utils.MergeBams as t_49_MergeLongbowPassedReclaimable_S2e { input: bams = t_17_FilterS2EReclaimableReads.passed_reads, prefix = SM + "_ccs_rejected_reclaimable_annotated_longbow_passed" }
    call Utils.MergeBams as t_50_MergeLongbowFailedReclaimable_S2e { input: bams = t_17_FilterS2EReclaimableReads.failed_reads, prefix = SM + "_ccs_rejected_reclaimable_annotated_longbow_failed" }

    # All Longbow Passed / Failed reads:
    call Utils.MergeBams as t_51_MergeAllLongbowPassedReads_S2e { input: bams = t_18_MergeLongbowS2EPassedReads.merged_bam, prefix = SM + "_all_longbow_passed" }
    call Utils.MergeBams as t_52_MergeAllLongbowFailedReads_S2e { input: bams = t_19_MergeLongbowS2EFailedReads.merged_bam, prefix = SM + "_all_longbow_failed" }

    # Merge CCS Reclaimed Array elements:
    call Utils.MergeBams as t_53_MergeCCSReclaimedArrayElements_S2e { input: bams = t_24_MergeS2ECcsReclaimedArrayElementSubshards.merged_bam, prefix = SM + "_ccs_reclaimed_array_elements"  }


    # Alias out the data we need to pass into stuff later:
    File ccs_corrected_reads = t_42_MergeCCSRqFilteredReads.merged_bam
    File ccs_corrected_reads_index = t_42_MergeCCSRqFilteredReads.merged_bai
    File ccs_rejected_reads = t_43_MergeCCSRqRejectedReads.merged_bam
    File ccs_rejected_reads_index = t_43_MergeCCSRqRejectedReads.merged_bai
    File annotated_ccs_reads = t_44_MergeAnnotatedCCSReads_S2e.merged_bam
    File annotated_ccs_reads_index = t_44_MergeAnnotatedCCSReads_S2e.merged_bai
    File longbow_passed_ccs_reads = t_45_MergeLongbowPassedCCSReads_S2e.merged_bam
    File longbow_passed_ccs_reads_index = t_45_MergeLongbowPassedCCSReads_S2e.merged_bai
    File longbow_failed_ccs_reads = t_46_MergeLongbowFailedCCSReads_S2e.merged_bam
    File longbow_failed_ccs_reads_index = t_46_MergeLongbowFailedCCSReads_S2e.merged_bai
    File ccs_reclaimable_reads = t_47_MergeCCSReclaimableReads_S2e.merged_bam
    File ccs_reclaimable_reads_index = t_47_MergeCCSReclaimableReads_S2e.merged_bai
    
    File annotated_ccs_reclaimable_reads = t_48_MergeCCSReclaimableAnnotatedReads_S2e.merged_bam
    
    File annotated_ccs_reclaimable_reads_index = t_48_MergeCCSReclaimableAnnotatedReads_S2e.merged_bai
    File ccs_reclaimed_reads = t_49_MergeLongbowPassedReclaimable_S2e.merged_bam
    File ccs_reclaimed_reads_index = t_49_MergeLongbowPassedReclaimable_S2e.merged_bai
    File longbow_failed_ccs_unreclaimable_reads = t_50_MergeLongbowFailedReclaimable_S2e.merged_bam
    File longbow_failed_ccs_unreclaimable_reads_index = t_50_MergeLongbowFailedReclaimable_S2e.merged_bai
    File longbow_passed_reads = t_51_MergeAllLongbowPassedReads_S2e.merged_bam
    File longbow_passed_reads_index = t_51_MergeAllLongbowPassedReads_S2e.merged_bai
    File longbow_failed_reads = t_52_MergeAllLongbowFailedReads_S2e.merged_bam
    File longbow_failed_reads_index = t_52_MergeAllLongbowFailedReads_S2e.merged_bai

    File ccs_reclaimed_array_elements = t_53_MergeCCSReclaimedArrayElements_S2e.merged_bam
    File ccs_reclaimed_array_elements_index = t_53_MergeCCSReclaimedArrayElements_S2e.merged_bai

    # Merge all CCS bams together for this Subread BAM:
    RuntimeAttr merge_extra_cpu_attrs = object {
        cpu_cores: 4
    }
    call Utils.MergeBams as t_54_MergeLongbowExtractedArrayElements { input: bams = t_39_ExtractCodingRegionsFromArrayElements.extracted_reads, prefix = SM + "_array_elements_longbow_extracted" }
    call Utils.MergeBams as t_55_MergeGenomeAlignedExtractedArrayElements { input: bams = t_41_RestoreAnnotationsToGenomeAlignedBam.output_bam, prefix = SM + "_array_elements_longbow_extracted_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }

    # We must discover the transcriptome for this new sample.
    # We only do this if we don't have SIRV data:
    if ( !is_SIRV_data ) {
        call StringTie2.Quantify as t_56_ST2_Quant {
            input:
                aligned_bam = t_55_MergeGenomeAlignedExtractedArrayElements.merged_bam,
                aligned_bai = t_55_MergeGenomeAlignedExtractedArrayElements.merged_bai,
                gtf = genome_annotation_gtf,
                keep_retained_introns = false,
                prefix = SM + "_StringTie2_Quantify",
        }

        call StringTie2.ExtractTranscriptSequences as t_57_ST2_ExtractTranscriptSequences  {
            input:
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_index,
                gtf = t_56_ST2_Quant.st_gtf,
                prefix = SM + "_StringTie2_ExtractTranscriptSequences",
        }

        call StringTie2.CompareTranscriptomes as t_58_ST2_CompareTranscriptomes {
            input:
                guide_gtf = genome_annotation_gtf,
                new_gtf = t_56_ST2_Quant.st_gtf,
                prefix = SM + "_StringTie2_CompareTranscriptome",
        }
    }

    # Set our transcriptome files:
    File transcriptome_reference_for_quant = if is_SIRV_data then transcriptome_ref_fasta else select_first([t_57_ST2_ExtractTranscriptSequences.transcripts_fa])
    File transcriptome_reference_index_for_quant = if is_SIRV_data then transcriptome_ref_fasta_index else select_first([t_57_ST2_ExtractTranscriptSequences.transcripts_fai])
    File transcriptome_reference_dict_for_quant = if is_SIRV_data then transcriptome_ref_fasta_dict else select_first([t_57_ST2_ExtractTranscriptSequences.transcripts_dict])

    # Now we have to align the array elements to the new transcriptome.
    scatter (extracted_array_elements in t_39_ExtractCodingRegionsFromArrayElements.extracted_reads) {
        # Align our array elements:
        call AR.Minimap2 as t_59_AlignArrayElementsToTranscriptome {
            input:
                reads      = [ extracted_array_elements ],
                ref_fasta  = transcriptome_reference_for_quant,
                map_preset = "asm20"
        }
        # We need to restore the annotations we created with the 10x tool to the aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_60_RestoreAnnotationsToTranscriptomeAlignedBam {
            input:
                annotated_bam_file = extracted_array_elements,
                aligned_bam_file = t_59_AlignArrayElementsToTranscriptome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        # To properly count our transcripts we must throw away the non-primary and unaligned reads:
        RuntimeAttr filterReadsAttrs = object {
            cpu_cores: 4,
            preemptible_tries: 0
        }
        call Utils.FilterReadsBySamFlags as t_61_RemoveUnmappedAndNonPrimaryReads {
            input:
                bam = t_60_RestoreAnnotationsToTranscriptomeAlignedBam.output_bam,
                sam_flags = "2308",
                prefix = SM + "_ArrayElements_Annotated_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }
        # Filter reads with no UMI tag:
        call Utils.FilterReadsWithTagValues as t_62_FilterReadsWithNoUMI {
            input:
                bam = t_61_RemoveUnmappedAndNonPrimaryReads.output_bam,
                tag = "ZU",
                value_to_remove = ".",
                prefix = SM + "_ArrayElements_Annotated_Aligned_PrimaryOnly_WithUMIs",
                runtime_attr_override = filterReadsAttrs
        }
        # Copy the contig to a tag.
        # By this point in the pipeline, array elements are aligned to a transcriptome, so this tag will
        # actually indicate the transcript to which each array element aligns.
        call TENX.CopyContigNameToReadTag as t_63_CopyContigNameToReadTag {
            input:
                aligned_bam_file = t_62_FilterReadsWithNoUMI.output_bam,
                prefix = SM + "_ArrayElements_Annotated_Aligned_PrimaryOnly_WithUMIs"
        }
    }

    # Now we merge together our TX-ome aligned stuff:
    call Utils.MergeBams as t_64_MergeTranscriptomeAlignedExtractedArrayElements { input: bams = t_60_RestoreAnnotationsToTranscriptomeAlignedBam.output_bam, prefix = SM + "_array_elements_longbow_extracted_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_65_MergePrimaryTranscriptomeAlignedArrayElements { input: bams = t_63_CopyContigNameToReadTag.output_bam, prefix = SM + "_array_elements_longbow_extracted_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    # Collect metrics on the subreads bam:
    RuntimeAttr subreads_sam_stats_runtime_attrs = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            ceil(3 * size(reads_bam, "GiB")),
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.8"
    }
    call AM.SamtoolsStats as t_66_CalcSamStatsOnInputBam {
        input:
            bam = reads_bam,
            runtime_attr_override = subreads_sam_stats_runtime_attrs
    }

    ##########
    # Quantify Transcripts:
    ##########

    call UMI_TOOLS.Run_Group as t_67_UMIToolsGroup {
        input:
            aligned_transcriptome_reads = t_65_MergePrimaryTranscriptomeAlignedArrayElements.merged_bam,
            aligned_transcriptome_reads_index = t_65_MergePrimaryTranscriptomeAlignedArrayElements.merged_bai,
            do_per_cell = !is_SIRV_data,
            prefix = "~{SM}_~{ID}_umi_tools_group"
    }

    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_68_CreateCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_67_UMIToolsGroup.output_bam,
            prefix = "~{SM}_~{ID}_gene_tx_expression_count_matrix"
    }

    # Only create the anndata objects if we're looking at real genomic data:
    if ( ! is_SIRV_data ) {
        call TX_POST.CreateCountMatrixAnndataFromTsv as t_69_CreateCountMatrixAnndataFromTsv {
            input:
                count_matrix_tsv = t_68_CreateCountMatrixFromAnnotatedBam.count_matrix,
                genome_annotation_gtf_file = select_first([t_56_ST2_Quant.st_gtf]),
                gencode_reference_gtf_file = genome_annotation_gtf,
                overlap_intervals = intervals_of_interest,
                overlap_interval_label = interval_overlap_name,
                prefix = "~{SM}_~{ID}_gene_tx_expression_count_matrix"
        }
    }

    ############################################################
    #               __  __      _        _
    #              |  \/  | ___| |_ _ __(_) ___ ___
    #              | |\/| |/ _ \ __| '__| |/ __/ __|
    #              | |  | |  __/ |_| |  | | (__\__ \
    #              |_|  |_|\___|\__|_|  |_|\___|___/
    #
    ############################################################

    String base_out_dir = outdir + "/" + DIR + "/" + t_01_WdlExecutionStartTimestamp.timestamp_string
    String metrics_out_dir = base_out_dir + "/metrics"

    # Aligned CCS Metrics:
    call RM.CalculateAndFinalizeReadMetrics as t_70_GenomeAlignedArrayElementMetrics {
        input:
            bam_file = t_55_MergeGenomeAlignedExtractedArrayElements.merged_bam,
            bam_index = t_55_MergeGenomeAlignedExtractedArrayElements.merged_bai,
            ref_dict = ref_fasta_dict,

            base_metrics_out_dir = metrics_out_dir + "/genome_aligned_array_element_metrics"
    }

    # Aligned Array Element Metrics:
    call RM.CalculateAndFinalizeAlternateReadMetrics as t_71_TranscriptomeAlignedArrayElementMetrics {
        input:
            bam_file = t_64_MergeTranscriptomeAlignedExtractedArrayElements.merged_bam,
            bam_index = t_64_MergeTranscriptomeAlignedExtractedArrayElements.merged_bai,
            ref_dict = transcriptome_reference_dict_for_quant,

            base_metrics_out_dir = metrics_out_dir + "/transcriptome_aligned_array_element_metrics"
    }

    # We should have exactly 1 CCS report for Sequel IIe reads:
    File final_ccs_report = t_10_FindCCSReport.ccs_report[0]

    ##########################################################################################
    #         ____                _                ____                       _
    #        / ___|_ __ ___  __ _| |_ ___         |  _ \ ___ _ __   ___  _ __| |_
    #       | |   | '__/ _ \/ _` | __/ _ \        | |_) / _ \ '_ \ / _ \| '__| __|
    #       | |___| | |  __/ (_| | ||  __/        |  _ <  __/ |_) | (_) | |  | |_
    #        \____|_|  \___|\__,_|\__\___|        |_| \_\___| .__/ \___/|_|   \__|
    #                                                       |_|
    ##########################################################################################

    RuntimeAttr create_report_runtime_attrs = object {
            preemptible_tries:  0
    }
    call JUPYTER.PB10xMasSeqSingleFlowcellReport as t_72_GenerateStaticReport {
        input:
            notebook_template                 = jupyter_template_static,

            sample_name                       = SM,

            subreads_stats                    = t_66_CalcSamStatsOnInputBam.raw_stats,
            ccs_reads_stats                   = t_70_GenomeAlignedArrayElementMetrics.sam_stats_raw_stats,
            array_elements_stats              = t_71_TranscriptomeAlignedArrayElementMetrics.sam_stats_raw_stats,
            ccs_report_file                   = final_ccs_report,

            raw_ccs_bam_file                  = ccs_corrected_reads,
            array_element_bam_file            = t_64_MergeTranscriptomeAlignedExtractedArrayElements.merged_bam,
            array_elements_genome_aligned     = t_55_MergeGenomeAlignedExtractedArrayElements.merged_bam,
            ccs_rejected_bam_file             = ccs_rejected_reads,

            annotated_bam_file                = annotated_ccs_reads,

            longbow_passed_reads_file         = longbow_passed_reads,
            longbow_failed_reads_file         = longbow_failed_reads,

            longbow_passed_ccs_reads          = longbow_passed_ccs_reads,
            longbow_failed_ccs_reads          = longbow_failed_ccs_reads,
            ccs_reclaimable_reads             = annotated_ccs_reclaimable_reads,
            ccs_reclaimed_reads               = ccs_reclaimed_reads,
            ccs_rejected_longbow_failed_reads = longbow_failed_ccs_unreclaimable_reads,
            raw_array_elements                = annotated_array_elements,
            ccs_reclaimed_array_elements      = ccs_reclaimed_array_elements,

            zmw_stats_json_gz                 = t_03_FindZmwStatsJsonGz.zmw_stats_json_gz,

#            ten_x_metrics_file                = t_28_TenxAnnotateArrayElements.stats,
            mas_seq_model                     = mas_seq_model,

            workflow_dot_file                 = workflow_dot_file,
            prefix                            = SM + "_MAS-seq_",

            runtime_attr_override             = create_report_runtime_attrs,
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

    String array_element_dir = base_out_dir + "/annotated_array_elements"
    String intermediate_reads_dir = base_out_dir + "/intermediate_reads"
    String quant_dir = base_out_dir + "/quant"
    String report_dir = base_out_dir + "/report"

    ##############################################################################################################
    # Finalize the final annotated, aligned array elements:
    call FF.FinalizeToDir as t_73_FinalizeQuantifiedArrayElements {
        input:
            files = [
                t_65_MergePrimaryTranscriptomeAlignedArrayElements.merged_bam,
                t_65_MergePrimaryTranscriptomeAlignedArrayElements.merged_bai,
#                PbIndexPrimaryTranscriptomeAlignedArrayElements.pbindex
            ],
            outdir = array_element_dir,
            keyfile = t_72_GenerateStaticReport.html_report
    }

    ##############################################################################################################
    # Finalize the discovered transcript:
    if ( !is_SIRV_data ) {
        call FF.FinalizeToDir as t_74_FinalizeDiscoveredTranscriptome {
            input:
                files = [
                    select_first([t_56_ST2_Quant.st_gtf]),
                    select_first([t_57_ST2_ExtractTranscriptSequences.transcripts_fa]),
                    select_first([t_57_ST2_ExtractTranscriptSequences.transcripts_fai]),
                    select_first([t_57_ST2_ExtractTranscriptSequences.transcripts_dict]),
                    select_first([t_58_ST2_CompareTranscriptomes.annotated_gtf]),
                    select_first([t_58_ST2_CompareTranscriptomes.loci]),
                    select_first([t_58_ST2_CompareTranscriptomes.stats]),
                    select_first([t_58_ST2_CompareTranscriptomes.tracking]),
                    select_first([t_58_ST2_CompareTranscriptomes.refmap]),
                    select_first([t_58_ST2_CompareTranscriptomes.tmap]),
                ],
                outdir = base_out_dir + "/discovered_transcriptome",
                keyfile = t_72_GenerateStaticReport.html_report
        }
    }
    ##############################################################################################################
    # Finalize the intermediate reads files (from raw CCS corrected reads through split array elements)
    call FF.FinalizeToDir as t_75_FinalizeArrayReads {
        input:
            files = [
                ccs_corrected_reads,
                ccs_corrected_reads_index,
                ccs_rejected_reads,
                ccs_rejected_reads_index,
                annotated_ccs_reads,
                annotated_ccs_reads_index,
                longbow_passed_ccs_reads,
                longbow_passed_ccs_reads_index,
                longbow_failed_ccs_reads,
                longbow_failed_ccs_reads_index,
                ccs_reclaimable_reads,
                ccs_reclaimable_reads_index,
                annotated_ccs_reclaimable_reads,
                annotated_ccs_reclaimable_reads_index,
                ccs_reclaimed_reads,
                ccs_reclaimed_reads_index,
                longbow_failed_ccs_unreclaimable_reads,
                longbow_failed_ccs_unreclaimable_reads_index,
                longbow_passed_reads,
                longbow_passed_reads_index,
                longbow_failed_reads,
                longbow_failed_reads_index
            ],
            outdir = intermediate_reads_dir + "/array_bams",
            keyfile = t_72_GenerateStaticReport.html_report
    }

    call FF.FinalizeToDir as t_76_FinalizeArrayElementReads {
        input:
            files = [
                annotated_array_elements,
                t_54_MergeLongbowExtractedArrayElements.merged_bam,
                t_54_MergeLongbowExtractedArrayElements.merged_bai,
                t_64_MergeTranscriptomeAlignedExtractedArrayElements.merged_bam,
                t_64_MergeTranscriptomeAlignedExtractedArrayElements.merged_bai,
                t_55_MergeGenomeAlignedExtractedArrayElements.merged_bam,
                t_55_MergeGenomeAlignedExtractedArrayElements.merged_bai,
            ],
            outdir = intermediate_reads_dir + "/array_element_bams",
            keyfile = t_72_GenerateStaticReport.html_report
    }

    ##############################################################################################################
    # Finalize Metrics:
    call FF.FinalizeToDir as t_77_FinalizeSamStatsOnInputBam {
        input:
            # an unfortunate hard-coded path here:
            outdir = metrics_out_dir + "/input_bam_stats",
            files = [
                t_66_CalcSamStatsOnInputBam.raw_stats,
                t_66_CalcSamStatsOnInputBam.summary_stats,
                t_66_CalcSamStatsOnInputBam.first_frag_qual,
                t_66_CalcSamStatsOnInputBam.last_frag_qual,
                t_66_CalcSamStatsOnInputBam.first_frag_gc_content,
                t_66_CalcSamStatsOnInputBam.last_frag_gc_content,
                t_66_CalcSamStatsOnInputBam.acgt_content_per_cycle,
                t_66_CalcSamStatsOnInputBam.insert_size,
                t_66_CalcSamStatsOnInputBam.read_length_dist,
                t_66_CalcSamStatsOnInputBam.indel_distribution,
                t_66_CalcSamStatsOnInputBam.indels_per_cycle,
                t_66_CalcSamStatsOnInputBam.coverage_distribution,
                t_66_CalcSamStatsOnInputBam.gc_depth
            ],
            keyfile = t_72_GenerateStaticReport.html_report
    }

    # Finalize all the 10x metrics here:
    # NOTE: We only run the 10x tool if we have real (non-SIRV) data, so we have to have this conditional here:
    if (! is_SIRV_data) {
        String tenXToolMetricsDir = metrics_out_dir + "/ten_x_tool_metrics"

        call FF.FinalizeToDir as t_78_FinalizeTenXRgStats {
            input:
                files = select_all([
                    starcode_seeds
               ]),
                outdir = tenXToolMetricsDir,
                keyfile = t_72_GenerateStaticReport.html_report
        }
    }

    call FF.FinalizeToDir as t_79_FinalizeCCSMetrics {
        input:
            files = [ final_ccs_report ],
            outdir = metrics_out_dir + "/ccs_metrics",
            keyfile = t_72_GenerateStaticReport.html_report
    }

    ##############################################################################################################
    # Finalize all the Quantification data:
    call FF.FinalizeToDir as t_80_FinalizeQuantResults {
        input:
            files = [
                t_67_UMIToolsGroup.output_bam,
                t_67_UMIToolsGroup.output_tsv,
                t_68_CreateCountMatrixFromAnnotatedBam.count_matrix
            ],
            outdir = quant_dir,
            keyfile = t_72_GenerateStaticReport.html_report
    }
    # Finalize our anndata objects if we have them:
    if ( ! is_SIRV_data ) {
        call FF.FinalizeToDir as t_81_FinalizeProcessedQuantResults {
            input:
                files = select_all([
                    t_69_CreateCountMatrixAnndataFromTsv.transcript_gene_count_anndata_h5ad,
                ]),
                outdir = quant_dir,
                keyfile = t_72_GenerateStaticReport.html_report
        }

        call FF.FinalizeToDir as t_82_FinalizeProcessedQuantResultsPickles {
            input:
                files = select_first([t_69_CreateCountMatrixAnndataFromTsv.pickles]),
                outdir = quant_dir,
                keyfile = t_72_GenerateStaticReport.html_report
        }
    }

    ##############################################################################################################
    # Finalize the report:
    call FF.FinalizeToDir as t_83_FinalizeStaticReport {
        input:
            files = [
                t_72_GenerateStaticReport.populated_notebook,
                t_72_GenerateStaticReport.html_report,
            ],
            outdir = report_dir,
            keyfile = t_72_GenerateStaticReport.html_report
    }

    call FF.FinalizeTarGzContents as t_84_FinalizeReportFigures {
        input:
            tar_gz_file = t_72_GenerateStaticReport.figures_tar_gz,
            outdir = report_dir,
            keyfile = t_72_GenerateStaticReport.html_report
    }

    call FF.FinalizeToDir as t_85_FinalizeReportPickles {
        input:
            files = t_72_GenerateStaticReport.pickles,
            outdir = report_dir,
            keyfile = t_72_GenerateStaticReport.html_report
    }

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good:
    call FF.WriteCompletionFile as t_86_WriteCompletionFile {
        input:
            outdir = base_out_dir + "/",
            keyfile = t_72_GenerateStaticReport.html_report
    }
}
