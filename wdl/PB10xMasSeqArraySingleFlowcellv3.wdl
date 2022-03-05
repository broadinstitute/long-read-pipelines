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

workflow PB10xMasSeqSingleFlowcellv3 {

    meta {
        description : "This workflow is designed to process data from the MASSeq v2 protocol and produce aligned reads that are ready for downstream analysis (e.g. transcript isoform identification).  It takes in a raw PacBio run folder location on GCS and produces a folder containing the aligned reads and other processed data."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String gcs_input_dir
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/PB10xMasSeqSingleFlowcellv3"

        File segments_fasta = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/cDNA_array_15x.unique_seqs_for_cartographer.fasta"
        File boundaries_file = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/bounds_file_for_extraction.txt"

        File head_adapter_fasta = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/10x_adapter.fasta"
        File tail_adapter_fasta = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/tso_adapter.fasta"
        File ten_x_cell_barcode_whitelist = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/737K-august-2016.txt"

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

        String starcode_extra_params = "--dist 2 --sphere"
#        String starcode_extra_params = "--dist 2 --cluster-ratio 10"

        String expanded_cbc_tag = "vx"

        File? illumina_barcoded_bam

        # Default here is 0 because ccs uncorrected reads all seem to have RQ = -1.
        # All pathologically long reads also have RQ = -1.
        # This way we preserve the vast majority of the data, even if it has low quality.
        # We can filter it out at later steps.
        Float min_read_quality = 0.0
        Int max_reclamation_length = 60000

        Boolean is_SIRV_data = false
        String mas_seq_model = "mas15v2"

        String? sample_name
    }

    parameter_meta {
        gcs_input_dir : "Input folder on GCS in which to search for BAM files to process."
        gcs_out_root_dir : "Root output GCS folder in which to place results of this workflow."

        segments_fasta : "FASTA file containing unique segments for which to search in the given BAM files.   These segments are used as delimiters in the reads.  Read splitting uses these delimiters and the boundaries file."
        boundaries_file : "Text file containing two comma-separated segment names from the segments_fasta on each line.  These entries define delimited sections to be extracted from the reads and treated as individual array elements."

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

    scatter (sharded_reads in t_09_ShardLongReads.unmapped_shards) {

        ## No more preemption on this sharding - takes too long otherwise.
        RuntimeAttr disable_preemption_runtime_attrs = object {
            preemptible_tries: 0
        }

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
                prefix = SM + "_ccs_reads_subshard",
                num_shards = 10,
        }
        scatter (s2e_ccs_longbow_passed_shard in t_23_ShardS2ECcsLongbowPassedReads.unmapped_shards) {
            # Segment CCS reads into array elements:
            call LONGBOW.Segment as t_24_SegmentS2ECcsReads {
                input:
                    annotated_reads = s2e_ccs_longbow_passed_shard,
                    prefix = SM + "_ccs_array_elements_subshard",
                    extra_args = "-b",
                    model = mas_seq_model
            }
        }

        # Merge Filtered CCS reads together:
        call Utils.MergeBams as t_25_MergeCCSArrayElementShards {
            input:
                bams = t_24_SegmentS2ECcsReads.segmented_bam,
                prefix = SM + "_ccs_array_elements_shard"
        }

        # Merge CCS Barcode Conf files:
        call Utils.MergeFiles as t_26_MergeCCSBarcodeConfShards {
            input:
                files_to_merge = t_24_SegmentS2ECcsReads.barcode_conf_file,
                merged_file_name = SM + "_ccs_array_element_barcode_confs_shard.txt"
        }

        #####################

        # Shard our CCS reads into smaller problems to do work on array elements:
        call PB.ShardLongReads as t_27_ShardS2ECcsReclaimedReads {
            input:
                unaligned_bam = t_18_FilterS2EReclaimableReads.passed_reads,
                unaligned_pbi = t_21_PbIndexS2ELongbowPassedReclaimedReads.pbindex,
                prefix = SM + "_ccs_reclaimed_reads_subshard",
                num_shards = 10,
        }
        scatter (s2e_ccs_reclaimed_shard in t_27_ShardS2ECcsReclaimedReads.unmapped_shards) {
             # Segment Reclaimed reads into array elements:
            call LONGBOW.Segment as t_28_SegmentS2ECcsReclaimedReads {
                input:
                    annotated_reads = s2e_ccs_reclaimed_shard,
                    prefix = SM + "_ccs_reclaimed_array_elements_subshard",
                    extra_args = "-b",
                    model = mas_seq_model
            }
        }

        # Merge Filtered CCS Reclaimed reads together:
        call Utils.MergeBams as t_29_MergeCCSReclaimedArrayElementShards {
            input:
                bams = t_28_SegmentS2ECcsReclaimedReads.segmented_bam,
                prefix = SM + "_ccs_reclaimed_array_elements_shard"
        }

        # Merge CCS Barcode Conf files:
        call Utils.MergeFiles as t_30_MergeCCSReclaimedBarcodeConfShards {
            input:
                files_to_merge = t_28_SegmentS2ECcsReclaimedReads.barcode_conf_file,
                merged_file_name = SM + "_ccs_reclaimed_array_element_barcode_confs_shard.txt"
        }

        ###############

        # Now align the array elements with their respective alignment presets:

        # Align CCS reads to the genome:
        call AR.Minimap2 as t_31_AlignCCSArrayElementsToGenome {
            input:
                reads      = [ t_25_MergeCCSArrayElementShards.merged_bam ],
                ref_fasta  = ref_fasta,
                map_preset = "splice:hq"
        }

        # Align Reclaimed reads to the genome:
        call AR.Minimap2 as t_32_AlignReclaimedArrayElementsToGenome {
            input:
                reads      = [ t_29_MergeCCSReclaimedArrayElementShards.merged_bam ],
                ref_fasta  = ref_fasta,
                map_preset = "splice"
        }

        ##############

        # Now restore the tags to the aligned bam files:
        call TENX.RestoreAnnotationstoAlignedBam as t_33_RestoreAnnotationsToGenomeAlignedCCSBam {
            input:
                annotated_bam_file = t_25_MergeCCSArrayElementShards.merged_bam,
                aligned_bam_file = t_31_AlignCCSArrayElementsToGenome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        call TENX.RestoreAnnotationstoAlignedBam as t_34_RestoreAnnotationsToGenomeAlignedReclaimedBam {
            input:
                annotated_bam_file = t_29_MergeCCSReclaimedArrayElementShards.merged_bam,
                aligned_bam_file = t_32_AlignReclaimedArrayElementsToGenome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }
    }

    # Merge the arrays:
    call Utils.MergeBams as t_35_MergeCCSLongbowPassedArrayReads {
        input:
            bams = t_17_FilterS2ECCSReads.passed_reads,
            prefix = SM + "_ccs_array_reads_longbow_passed"
    }
    call Utils.MergeBams as t_36_MergeCCSLongbowFailedArrayReads {
        input:
            bams = t_17_FilterS2ECCSReads.failed_reads,
            prefix = SM + "_ccs_array_reads_longbow_failed"
    }

    call Utils.MergeBams as t_37_MergeCCSReclaimedArrayReads {
        input:
            bams = t_18_FilterS2EReclaimableReads.passed_reads,
            prefix = SM + "_ccs_reclaimed_array_reads_longbow_passed"
    }
    call Utils.MergeBams as t_38_MergeCCSUnreclaimableArrayReads {
        input:
            bams = t_18_FilterS2EReclaimableReads.failed_reads,
            prefix = SM + "_ccs_unreclaimable_array_reads_longbow_failed"
    }

    # Merge Filtered CCS reads together:
    call Utils.MergeBams as t_39_MergeCCSArrayElements {
        input:
            bams = t_25_MergeCCSArrayElementShards.merged_bam,
            prefix = SM + "_ccs_array_elements"
    }

    # Merge Filtered CCS Reclaimed reads together:
    call Utils.MergeBams as t_40_MergeCCSReclaimedArrayElements {
        input:
            bams = t_29_MergeCCSReclaimedArrayElementShards.merged_bam,
            prefix = SM + "_ccs_reclaimed_array_elements"
    }

    # Merge Aligned CCS reads together:
    call Utils.MergeBams as t_41_MergeAlignedCCSArrayElements {
        input:
            bams = t_33_RestoreAnnotationsToGenomeAlignedCCSBam.output_bam,
            prefix = SM + "_ccs_array_elements_aligned"
    }

    # Merge Aligned CCS Reclaimed reads together:
    call Utils.MergeBams as t_42_MergeAlignedCCSReclaimedArrayElements {
        input:
            bams = t_34_RestoreAnnotationsToGenomeAlignedReclaimedBam.output_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_aligned"
    }

    # Merge CCS Barcode Conf files:
    call Utils.MergeFiles as t_43_MergeAllCCSBarcodeConfShards {
        input:
            files_to_merge = t_26_MergeCCSBarcodeConfShards.merged_file,
            merged_file_name = SM + "_ccs_array_element_barcode_confs.txt"
    }

    # Merge CCS Barcode Conf files:
    call Utils.MergeFiles as t_44_MergeAllCCSReclaimedBarcodeConfShards {
        input:
            files_to_merge = t_30_MergeCCSReclaimedBarcodeConfShards.merged_file,
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
    call Utils.FilterMasSeqReadsWithGatk as t_45_AlignmentFilterForCcsArrayElements {
        input:
            bam_file = t_41_MergeAlignedCCSArrayElements.merged_bam,
            bam_index = t_41_MergeAlignedCCSArrayElements.merged_bai,
            prefix = SM + "_CCS_ArrayElements_Annotated_Aligned_PrimaryOnly",
            runtime_attr_override = filterReadsAttrs
    }

    call Utils.FilterMasSeqReadsWithGatk as t_46_AlignmentFilterForReclaimedArrayElements {
        input:
            bam_file = t_42_MergeAlignedCCSReclaimedArrayElements.merged_bam,
            bam_index = t_42_MergeAlignedCCSReclaimedArrayElements.merged_bai,
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
    call Utils.MergeBams as t_47_MergeAllAlignedAndFilteredArrayElements {
        input:
            bams = [t_45_AlignmentFilterForCcsArrayElements.bam, t_46_AlignmentFilterForReclaimedArrayElements.bam],
            prefix = SM + "_all_array_elements_aligned"
    }

    call StringTie2.Quantify as t_48_ST2_Quant {
        input:
            aligned_bam = t_47_MergeAllAlignedAndFilteredArrayElements.merged_bam,
            aligned_bai = t_47_MergeAllAlignedAndFilteredArrayElements.merged_bai,
            gtf = genome_annotation_gtf,
            keep_retained_introns = false,
            prefix = SM + "_StringTie2_Quantify",
    }

    call StringTie2.ExtractTranscriptSequences as t_49_ST2_ExtractTranscriptSequences  {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_index,
            gtf = t_48_ST2_Quant.st_gtf,
            prefix = SM + "_StringTie2_ExtractTranscriptSequences",
    }

    call StringTie2.CompareTranscriptomes as t_50_ST2_CompareTranscriptomes {
        input:
            guide_gtf = genome_annotation_gtf,
            new_gtf = t_48_ST2_Quant.st_gtf,
            prefix = SM + "_StringTie2_CompareTranscriptome",
    }

    ##########################################################################################################################
    ##########################################################################################################################

    # Now we pad our barcodes and correct them:
    Int num_array_element_shards = 50

    # CCS
    call Utils.ShardReads as t_51_ShardS2ECcsArrayElements {
        input:
            bam = t_45_AlignmentFilterForCcsArrayElements.bam,
            bam_index = t_45_AlignmentFilterForCcsArrayElements.bai,
            prefix = SM + "_ccs_array_elements_subshard",
            num_shards = num_array_element_shards,
    }
    scatter (ccsi in range(length(t_51_ShardS2ECcsArrayElements.shards))) {
        File ccs_array_element_shard = t_51_ShardS2ECcsArrayElements.shards[ccsi]

        call Utils.IndexBam as t_52_IndexCcsArrayElementShard {
            input:
                bam = ccs_array_element_shard
        }

        # Now that we've annotated the reads, we can pad the UMIs by a couple of bases to aid in the deduping:
        call LONGBOW.Pad as t_53_LongbowPadCCSArrayElementUMIs {
            input:
                reads = ccs_array_element_shard,
                model = mas_seq_model,
                tag_to_expand = "ZU",
                padding = 2,
                prefix = SM + "_ccs_array_elements_aligned_annotated_umi_padded_shard_" + ccsi
        }

        call LONGBOW.Pad as t_54_LongbowPadCCSArrayElementCBCs {
            input:
                reads = t_53_LongbowPadCCSArrayElementUMIs.padded_tag_bam,
                model = mas_seq_model,
                tag_to_expand = "CR",
                new_tag_dest = expanded_cbc_tag,
                padding = 3,
                prefix = SM + "_ccs_array_elements_aligned_annotated_cbc_padded_shard_" + ccsi
        }

        # Now we should correct our barcodes based on the whitelist:
        call LONGBOW.Correct as t_55_LongbowCorrectCCSCorrectedArrayElementCBCs {
            input:
                reads = t_54_LongbowPadCCSArrayElementCBCs.padded_tag_bam,
                barcode_allow_list = ten_x_cell_barcode_whitelist,
                model = mas_seq_model,
                ccs_lev_dist_threshold = 2,
                clr_lev_dist_threshold = 3,
                prefix = SM + "_ccs_array_elements_aligned_annotated_padded_cbc_corrected_shard_" + ccsi,
                raw_barcode_tag = expanded_cbc_tag,
                corrected_barcode_tag = "CB",
        }
    }
    # Merge Aligned CCS reads together:
    call Utils.MergeBams as t_56_MergeLongbowPaddedCBCCorrectedCCSArrayElements {
        input:
            bams = t_55_LongbowCorrectCCSCorrectedArrayElementCBCs.corrected_barcodes_bam,
            prefix = SM + "_ccs_array_elements_aligned_annotated_padded_CBC_corrected"
    }
    call Utils.MergeBams as t_57_MergeLongbowPaddedCBCUncorrectableCCSArrayElements {
        input:
            bams = t_55_LongbowCorrectCCSCorrectedArrayElementCBCs.uncorrected_barcodes_bam,
            prefix = SM + "_ccs_array_elements_aligned_annotated_padded_CBC_uncorrectable"
    }

    # RECLAIMED
    call Utils.ShardReads as t_58_ShardS2ECcsReclaimedArrayElements {
        input:
            bam = t_46_AlignmentFilterForReclaimedArrayElements.bam,
            bam_index = t_46_AlignmentFilterForReclaimedArrayElements.bai,
            prefix = SM + "_ccs_reclaimed_array_elements_subshard",
            num_shards = num_array_element_shards,
    }
    scatter (cri in range(length(t_58_ShardS2ECcsReclaimedArrayElements.shards))) {
        File ccs_reclaimed_array_element_shard = t_58_ShardS2ECcsReclaimedArrayElements.shards[cri]

        call Utils.IndexBam as t_59_IndexCcsReclaimedArrayElementShard {
            input:
                bam = ccs_reclaimed_array_element_shard
        }

        # Now that we've annotated the reads, we can pad the CBCs and UMIs by a couple of bases to aid in the deduping:
        call LONGBOW.Pad as t_60_LongbowPadCCSReclaimedArrayElementUMIs {
            input:
                reads = ccs_reclaimed_array_element_shard,
                model = mas_seq_model,
                tag_to_expand = "ZU",
                new_tag_dest = expanded_cbc_tag,
                padding = 2,
                prefix = SM + "_ccs_reclaimed_array_elements_aligned_annotated_umi_padded_shard_" + cri
        }

        call LONGBOW.Pad as t_61_LongbowPadCCSReclaimedArrayElementCBCs {
            input:
                reads = t_60_LongbowPadCCSReclaimedArrayElementUMIs.padded_tag_bam,
                model = mas_seq_model,
                tag_to_expand = "CR",
                padding = 3,
                prefix = SM + "_ccs_reclaimed_array_elements_aligned_annotated_cbc_padded_shard_" + cri
        }

        # Now we should correct our barcodes based on the whitelist:
        call LONGBOW.Correct as t_62_LongbowCorrectCCSReclaimedArrayElementCBCs {
            input:
                reads = t_61_LongbowPadCCSReclaimedArrayElementCBCs.padded_tag_bam,
                barcode_allow_list = ten_x_cell_barcode_whitelist,
                model = mas_seq_model,
                ccs_lev_dist_threshold = 2,
                clr_lev_dist_threshold = 3,
                prefix = SM + "_ccs_reclaimed_array_elements_aligned_annotated_padded_cbc_corrected_shard_" + cri,
                raw_barcode_tag = expanded_cbc_tag,
                corrected_barcode_tag = "CB",
        }
    }
    # Merge Aligned CCS Reclaimed reads together:
    call Utils.MergeBams as t_63_MergeLongbowPaddedCBCCorrectedCCSReclaimedArrayElements {
        input:
            bams = t_62_LongbowCorrectCCSReclaimedArrayElementCBCs.corrected_barcodes_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_aligned_annotated_padded_CBC_corrected"
    }
    call Utils.MergeBams as t_64_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElements {
        input:
            bams = t_62_LongbowCorrectCCSReclaimedArrayElementCBCs.uncorrected_barcodes_bam,
            prefix = SM + "_ccs_array_elements_aligned_annotated_padded_CBC_uncorrectable"
    }

    #####################
    # Now we correct the barcodes:
    #     ___  _     ____
    #    / _ \| |   |  _ \
    #   | | | | |   | | | |
    #   | |_| | |___| |_| |
    #    \___/|_____|____/

    # We need to merge ALL the barcode information together before correcting either the
    # CCS or the Reclaimed reads.

#    call Utils.MergeFiles as t_59_MergeCbcConfScoreTsvsForStarcode {
#        input:
#            files_to_merge = [t_43_MergeAllCCSBarcodeConfShards.merged_file, t_44_MergeAllCCSReclaimedBarcodeConfShards.merged_file],
#            merged_file_name = "all_array_element_raw_starcode_counts.txt"
#    }
#
#    # If we have our ilmn barcode file, we need to process it here:
#    if (defined(illumina_barcoded_bam)) {
#        call TENX.ExtractIlmnBarcodeConfScores as t_60_ExtractIlmnBarcodeConfScores {
#            input:
#                bam_file = select_first([illumina_barcoded_bam]),
#                prefix = SM,
#                runtime_attr_override = fast_network_attrs
#        }
#
#        # Concatenate the TSV files with the barcode scores that we just created:
#        call Utils.MergeFiles as t_61_GetMasterCbcConfScoreTsvForStarcode {
#            input:
#                files_to_merge = [t_59_MergeCbcConfScoreTsvsForStarcode.merged_file, t_60_ExtractIlmnBarcodeConfScores.conf_score_tsv],
#                merged_file_name = "combined_mas-seq_and_ilmn_raw_starcode_counts.txt"
#        }
#    }
#    File starcode_seeds = if (defined(illumina_barcoded_bam)) then select_first([t_61_GetMasterCbcConfScoreTsvForStarcode.merged_file]) else t_59_MergeCbcConfScoreTsvsForStarcode.merged_file
#
#    # We have to consolidate our seeds into unique entries for starcode not to crash and burn:
#    call TX_POST.MergeBarcodeCounts as t_62_ConsolidateBarcodeCountsForStarcode {
#        input:
#            barcode_count_tsv = starcode_seeds,
#            prefix = SM + "_barcode_counts_for_starcode"
#    }
#
#    # Now we can correct our barcodes:
#    call TENX.CorrectBarcodesWithStarcodeSeedCounts as t_63_CorrectCCSBarcodesWithStarcodeSeedCountsSharded {
#        input:
#            bam_file = t_54_MergeLongbowPaddedCCSArrayElements.merged_bam,
#            starcode_seeds_tsv = t_62_ConsolidateBarcodeCountsForStarcode.merged_counts,
#            whitelist_10x = ten_x_cell_barcode_whitelist,
#            extra_parameters = starcode_extra_params,
#            prefix = SM + "_annotated_ccs_array_elements_starcode_padded"
#    }
#    call TENX.CorrectBarcodesWithStarcodeSeedCounts as t_64_CorrectReclaimedBarcodesWithStarcodeSeedCountsSharded {
#        input:
#            bam_file = t_58_MergeLongbowPaddedCCSReclaimedArrayElements.merged_bam,
#            starcode_seeds_tsv = t_62_ConsolidateBarcodeCountsForStarcode.merged_counts,
#            whitelist_10x = ten_x_cell_barcode_whitelist,
#            extra_parameters = starcode_extra_params,
#            prefix = SM + "_annotated_ccs_reclaimed_array_elements_starcode_padded"
#    }
#
#    # Merge Aligned CCS and Reclaimed reads together:
#    call Utils.MergeBams as t_65_MergeAllAnnotatedArrayElements {
#        input:
#            bams = [t_63_CorrectCCSBarcodesWithStarcodeSeedCountsSharded.output_bam, t_64_CorrectReclaimedBarcodesWithStarcodeSeedCountsSharded.output_bam],
#            prefix = SM + "_all_starcode_annotated_array_elements_padded"
#    }


    #################
    # Here we restore the original read names to the bam because we're hashing them with Longbow.segment:

    # Merge Aligned CCS and Reclaimed reads together:
    call Utils.MergeBams as t_65_MergeAllAnnotatedArrayElements {
        input:
            bams = [t_56_MergeLongbowPaddedCBCCorrectedCCSArrayElements.merged_bam, t_63_MergeLongbowPaddedCBCCorrectedCCSReclaimedArrayElements.merged_bam],
            prefix = SM + "_all_cbc_annotated_array_elements_padded"
    }

    call TX_PRE.RestoreOriginalReadNames as t_66_RestoreOriginalReadNames {
        input:
            bam = t_65_MergeAllAnnotatedArrayElements.merged_bam,
            prefix =  SM + "_all_cbc_annotated_array_elements_padded_original_names"
    }

    #################
    # Now we have to split the reads again, process them into gff files, run gffcompare and then aggregate the results in a graph

    # We can actually compare the references without needing to scatter:
    call TX_PRE.GffCompare as t_67_GffCompareStringtie2toGencode {
        input:
            gff_ref = t_48_ST2_Quant.st_gtf,
            gff_query = genome_annotation_gtf,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
    }
    call TX_PRE.GffCompare as t_68_GffCompareGencodetoStringtie2 {
        input:
            gff_ref = genome_annotation_gtf,
            gff_query = t_48_ST2_Quant.st_gtf,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
    }

    # Split by contig:
    call TX_PRE.SplitBamByContig as t_69_SplitArrayElementsByContig {
        input:
            bam = t_66_RestoreOriginalReadNames.bam_out,
            prefix = SM + "_all_cbc_annotated_array_elements_padded_original_names"
    }

    # For each contig:
    scatter (i in range(length(t_69_SplitArrayElementsByContig.contig_bams))) {

        File contig_bam = t_69_SplitArrayElementsByContig.contig_bams[i]
        String contig_name = t_69_SplitArrayElementsByContig.contig_names[i]

        # Create a GFF file:
        call TX_PRE.ConvertSplicedBamToGff as t_70_ConvertSplicedBamToGff {
            input:
                bam = contig_bam
        }

        # Compare GFF files:
        call TX_PRE.GffCompare as t_71_GffCompareStringtie2toMasSeqReads {
            input:
                gff_ref = t_48_ST2_Quant.st_gtf,
                gff_query = t_70_ConvertSplicedBamToGff.gff,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
        }

        call TX_PRE.GffCompare as t_72_GffCompareGencodetoMasSeqReads {
            input:
                gff_ref = genome_annotation_gtf,
                gff_query = t_70_ConvertSplicedBamToGff.gff,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
        }

        # Create the comparison graph and tsv files:
        call TX_POST.QuantifyGffComparison as t_73_QuantifyGffComparison {
            input:
                genome_gtf = genome_annotation_gtf,
                st2_gencode_refmap = t_67_GffCompareStringtie2toGencode.refmap,
                st2_gencode_tmap = t_67_GffCompareStringtie2toGencode.tmap,
                st2_read_refmap = t_71_GffCompareStringtie2toMasSeqReads.refmap,
                st2_read_tmap = t_71_GffCompareStringtie2toMasSeqReads.tmap,
                gencode_st2_refmap = t_68_GffCompareGencodetoStringtie2.refmap,
                gencode_st2_tmap = t_68_GffCompareGencodetoStringtie2.tmap,
                gencode_read_refmap = t_72_GffCompareGencodetoMasSeqReads.refmap,
                gencode_read_tmap = t_72_GffCompareGencodetoMasSeqReads.tmap,
                prefix = SM + "_all_cbc_annotated_array_elements_padded_" + contig_name
        }
    }

    # Merge our tx equivalance classes assignments and eq classes:
    call TX_POST.CombineEqClassFiles as t_74_CombineEqClassFiles {
        input:
            gene_eq_class_definitions = t_73_QuantifyGffComparison.gene_eq_class_labels_file,
            gene_assignment_files = t_73_QuantifyGffComparison.gene_assignments_file,
            equivalence_class_definitions = t_73_QuantifyGffComparison.tx_equivalence_class_labels_file,
            equivalence_classes = t_73_QuantifyGffComparison.tx_equivalence_class_file,
            prefix = SM + "_all_cbc_annotated_array_elements_padded"
    }

    ############################################################
    # Quantify Transcripts:
    ##########

    # TODO:
    # Slot in Vic's fix for UMI assignments / dedupe here.

    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_75_CreateCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_66_RestoreOriginalReadNames.bam_out,
            tx_equivalence_class_assignments = t_74_CombineEqClassFiles.combined_tx_eq_class_assignments,
            umi_tag = "ZU",
            prefix = SM + "_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_76_CreateCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_75_CreateCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = t_48_ST2_Quant.st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = t_74_CombineEqClassFiles.combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = t_74_CombineEqClassFiles.combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = t_74_CombineEqClassFiles.combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = t_74_CombineEqClassFiles.combined_gene_eq_class_defs,
            prefix = SM + "_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: 64}
    }


    ############################################################
    #               __  __      _        _
    #              |  \/  | ___| |_ _ __(_) ___ ___
    #              | |\/| |/ _ \ __| '__| |/ __/ __|
    #              | |  | |  __/ |_| |  | | (__\__ \
    #              |_|  |_|\___|\__|_|  |_|\___|___/
    #
    ############################################################

    call LONGBOW.AggregateCorrectLogStats as t_77_AggregateLongbowCorrectStats {
        input:
            longbow_correct_log_files = flatten([t_62_LongbowCorrectCCSReclaimedArrayElementCBCs.log, t_55_LongbowCorrectCCSCorrectedArrayElementCBCs.log]),
            out_name = SM + "_longbow_correct_stats.txt"
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

    File keyfile = t_76_CreateCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad

    String base_out_dir = outdir + "/" + DIR + "/" + t_01_WdlExecutionStartTimestamp.timestamp_string
    String stats_out_dir = base_out_dir + "/stats"
    String array_element_dir = base_out_dir + "/annotated_array_elements"
    String intermediate_reads_dir = base_out_dir + "/intermediate_reads"

    String meta_files_dir = base_out_dir + "/meta_files"

    String intermediate_array_reads_dir = intermediate_reads_dir + "/array_reads"
    String intermediate_array_elements_dir = intermediate_reads_dir + "/array_elements"

    String quant_dir = base_out_dir + "/quant"

    ##############################################################################################################
    # Finalize gene / tx assignments:
    call FF.FinalizeToDir as t_78_FinalizeTxAndGeneAssignments {
        input:
            files = [
                t_74_CombineEqClassFiles.combined_gene_eq_class_defs,
                t_74_CombineEqClassFiles.combined_gene_eq_class_assignments,
                t_74_CombineEqClassFiles.combined_tx_eq_class_defs,
                t_74_CombineEqClassFiles.combined_tx_eq_class_assignments,
                t_75_CreateCountMatrixFromAnnotatedBam.count_matrix,
                t_76_CreateCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad,
            ],
            outdir = quant_dir,
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_79_FinalizeRawQuantPickles {
        input:
            files = t_76_CreateCountMatrixAnndataFromEqClasses.pickles,
            outdir = quant_dir,
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_80_FinalizeRefAndSt2Comparisons {
        input:
            files = [
                t_67_GffCompareStringtie2toGencode.refmap,
                t_67_GffCompareStringtie2toGencode.tmap,
                t_67_GffCompareStringtie2toGencode.tracking,
                t_67_GffCompareStringtie2toGencode.loci,
                t_67_GffCompareStringtie2toGencode.annotated_gtf,
                t_67_GffCompareStringtie2toGencode.stats,
                t_67_GffCompareStringtie2toGencode.log,

                t_68_GffCompareGencodetoStringtie2.refmap,
                t_68_GffCompareGencodetoStringtie2.tmap,
                t_68_GffCompareGencodetoStringtie2.tracking,
                t_68_GffCompareGencodetoStringtie2.loci,
                t_68_GffCompareGencodetoStringtie2.annotated_gtf,
                t_68_GffCompareGencodetoStringtie2.stats,
                t_68_GffCompareGencodetoStringtie2.log,
            ],
            outdir = quant_dir + "/gencode_and_stringtie2",
            keyfile = keyfile
    }

    # Finalize gene / tx assignment by contig:
    # NOTE: According to the scatter/gather documentation in the WDL spec, this will work correctly
    #       (https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#scatter--gather)
    scatter (i in range(length(t_69_SplitArrayElementsByContig.contig_bams))) {
        String contig = t_69_SplitArrayElementsByContig.contig_names[i]

        call FF.FinalizeToDir as t_81_FinalizeTxAndGeneAssignmentsByContig {
            input:
                files = [
                    t_71_GffCompareStringtie2toMasSeqReads.refmap[i],
                    t_71_GffCompareStringtie2toMasSeqReads.tmap[i],
                    t_71_GffCompareStringtie2toMasSeqReads.tracking[i],
                    t_71_GffCompareStringtie2toMasSeqReads.loci[i],
                    t_71_GffCompareStringtie2toMasSeqReads.annotated_gtf[i],
                    t_71_GffCompareStringtie2toMasSeqReads.stats[i],
                    t_71_GffCompareStringtie2toMasSeqReads.log[i],

                    t_72_GffCompareGencodetoMasSeqReads.refmap[i],
                    t_72_GffCompareGencodetoMasSeqReads.tmap[i],
                    t_72_GffCompareGencodetoMasSeqReads.tracking[i],
                    t_72_GffCompareGencodetoMasSeqReads.loci[i],
                    t_72_GffCompareGencodetoMasSeqReads.annotated_gtf[i],
                    t_72_GffCompareGencodetoMasSeqReads.stats[i],
                    t_72_GffCompareGencodetoMasSeqReads.log[i],

                    t_73_QuantifyGffComparison.gene_assignments_file[i],
                    t_73_QuantifyGffComparison.tx_equivalence_class_labels_file[i],
                    t_73_QuantifyGffComparison.tx_equivalence_class_file[i],
                    t_73_QuantifyGffComparison.graph_gpickle[i],
                ],
                outdir = quant_dir + "/by_contig/" + contig,
                keyfile = keyfile
        }
    }

    ##############################################################################################################
    # Finalize annotated, aligned array elements:
    call FF.FinalizeToDir as t_82_FinalizeIntermediateCBCAnnotatedArrayElements {
        input:
            files = [
                t_56_MergeLongbowPaddedCBCCorrectedCCSArrayElements.merged_bam,
                t_57_MergeLongbowPaddedCBCUncorrectableCCSArrayElements.merged_bam,
                t_63_MergeLongbowPaddedCBCCorrectedCCSReclaimedArrayElements.merged_bam,
                t_64_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElements.merged_bam,
                t_65_MergeAllAnnotatedArrayElements.merged_bam,
                t_65_MergeAllAnnotatedArrayElements.merged_bai,
            ],
            outdir = intermediate_array_elements_dir,
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_83_FinalizeCBCAnnotatedArrayElements {
        input:
            files = [
                t_66_RestoreOriginalReadNames.bam_out,
            ],
            outdir = array_element_dir,
            keyfile = keyfile
    }

    ##############################################################################################################
    # Finalize meta files:
    call FF.FinalizeToDir as t_84_FinalizeMeta {
        input:
            files = [
#                starcode_seeds,
#                t_62_ConsolidateBarcodeCountsForStarcode.merged_counts,
                ten_x_cell_barcode_whitelist,
                t_43_MergeAllCCSBarcodeConfShards.merged_file,
                t_44_MergeAllCCSReclaimedBarcodeConfShards.merged_file
            ],
            outdir = meta_files_dir,
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_85_FinalizeCCSCBCcorrectionLogsToMeta {
        input:
            files = t_55_LongbowCorrectCCSCorrectedArrayElementCBCs.log,
            outdir = meta_files_dir + "/" + "ccs_cbc_correction_logs",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_86_FinalizeCCSRejectedCBCcorrectionLogsToMeta {
        input:
            files = t_62_LongbowCorrectCCSReclaimedArrayElementCBCs.log,
            outdir = meta_files_dir + "/" + "ccs_rejected_cbc_correction_logs",
            keyfile = keyfile
    }

#    if (defined(illumina_barcoded_bam)) {
#        call FF.FinalizeToDir as t_84_FinalizeMetaIlmnBarcodeConfs {
#            input:
#                files = select_all([
#                    t_60_ExtractIlmnBarcodeConfScores.conf_score_tsv
#                ]),
#                outdir = meta_files_dir,
#                keyfile = keyfile
#        }
#    }

    ##############################################################################################################
    # Finalize the discovered transcriptome:
    if ( !is_SIRV_data ) {
        call FF.FinalizeToDir as t_87_FinalizeDiscoveredTranscriptome {
            input:
                files = [
                    t_48_ST2_Quant.st_gtf,
                    t_49_ST2_ExtractTranscriptSequences.transcripts_fa,
                    t_49_ST2_ExtractTranscriptSequences.transcripts_fai,
                    t_49_ST2_ExtractTranscriptSequences.transcripts_dict,
                    t_50_ST2_CompareTranscriptomes.annotated_gtf,
                    t_50_ST2_CompareTranscriptomes.loci,
                    t_50_ST2_CompareTranscriptomes.stats,
                    t_50_ST2_CompareTranscriptomes.tracking,
                    t_50_ST2_CompareTranscriptomes.refmap,
                    t_50_ST2_CompareTranscriptomes.tmap,
                ],
                outdir = base_out_dir + "/discovered_transcriptome",
                keyfile = keyfile
        }
    }
    ##############################################################################################################
    # Finalize the intermediate reads files (from raw CCS corrected reads through split array elements)
    call FF.FinalizeToDir as t_88_FinalizeArrayReads {
        input:
            files = [
                t_35_MergeCCSLongbowPassedArrayReads.merged_bam,
                t_36_MergeCCSLongbowFailedArrayReads.merged_bam,
                t_37_MergeCCSReclaimedArrayReads.merged_bam,
                t_38_MergeCCSUnreclaimableArrayReads.merged_bam,
            ],
            outdir = intermediate_reads_dir + "/array_bams",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_89_FinalizeCCSMetrics {
        input:
            files = [ t_11_FindCCSReport.ccs_report[0], t_77_AggregateLongbowCorrectStats.stats ],
            outdir = stats_out_dir + "/ccs_metrics",
            keyfile = keyfile
    }

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good:
    call FF.WriteCompletionFile as t_90_WriteCompletionFile {
        input:
            outdir = base_out_dir + "/",
            keyfile = keyfile
    }
}
