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

workflow PB10xMasSeqSingleFlowcellv5 {

    meta {
        description : "This workflow is designed to process data from the MASSeq v2 protocol and produce aligned reads that are ready for downstream analysis (e.g. transcript isoform identification).  It takes in a raw PacBio run folder location on GCS and produces a folder containing the aligned reads and other processed data."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String gcs_input_dir
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/PB10xMasSeqSingleFlowcellv5"

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
        String mas_seq_model = "mas_15+sc_10x5p"

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
    call Utils.GetCurrentTimestampString as t_001_WdlExecutionStartTimestamp { input: }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams as t_002_FindBams { input: gcs_input_dir = gcs_input_dir }
    call PB.FindZmwStatsJsonGz as t_003_FindZmwStatsJsonGz { input: gcs_input_dir = gcs_input_dir }

    # Check here if we found ccs bams or subread bams:
    Boolean use_subreads = t_002_FindBams.has_subreads
    Array[String] top_level_bam_files = if use_subreads then t_002_FindBams.subread_bams else t_002_FindBams.ccs_bams

    # Make sure we have **EXACTLY** one bam file to run on:
    if (length(top_level_bam_files) != 1) {
        call Utils.FailWithWarning as t_004_WARN1 { input: warning = "Error: Multiple BAM files found.  Cannot continue!" }
    }

    if (use_subreads) {
        call Utils.FailWithWarning as t_005_WARN2 { input: warning = "Error: This workflow now only supports data from the Sequel IIe." }
    }

    # Alias our bam file so we can work with it easier:
    File reads_bam = top_level_bam_files[0]

    call PB.GetPbReadGroupInfo as t_006_GetReadGroupInfo { input: gcs_bam_path = reads_bam }
    call PB.GetRunInfo as t_007_GetRunInfo { input: subread_bam = reads_bam }

    String SM  = select_first([sample_name, t_007_GetRunInfo.run_info["SM"]])
    String PL  = "PACBIO"
    String PU  = t_007_GetRunInfo.run_info["PU"]
    String DT  = t_007_GetRunInfo.run_info["DT"]
    String ID  = PU
    String DS  = t_007_GetRunInfo.run_info["DS"]
    String DIR = SM + "." + ID

    String RG_subreads  = "@RG\\tID:~{ID}.subreads\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
    String RG_consensus = "@RG\\tID:~{ID}.consensus\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
    String RG_array_elements = "@RG\\tID:~{ID}.array_elements\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

    # Check to see if we need to annotate our reads:
    call LONGBOW.CheckForAnnotatedArrayReads as t_008_CheckForAnnotatedReads {
        input:
            bam = reads_bam
    }

    File read_pbi = sub(reads_bam, ".bam$", ".bam.pbi")
    call PB.ShardLongReads as t_009_ShardLongReads {
        input:
            unaligned_bam = reads_bam,
            unaligned_pbi = read_pbi,
            prefix = SM + "_shard",
            num_shards = 100,
    }

    ## No more preemption on this sharding - takes too long otherwise.
    RuntimeAttr disable_preemption_runtime_attrs = object {
        preemptible_tries: 0
    }

    Array[String] tags_to_preserve =  [ "CB", "JB", "JC", "JD", "JF", "JX", "RC", "RG", "SG", "XA", "XB", "XC", "XF", "XM", "XN", "XQ", "XU", "YC", "YG", "YK", "YN", "YP", "YQ", "YS", "YV", "ZS", "ZU", "ec", "fn", "ic", "im", "is", "it", "np", "pz", "rn", "rq", "sn", "we", "ws", "zm" ]

    scatter (main_shard_index in range(length(t_009_ShardLongReads.unmapped_shards))) {
        File sharded_reads = t_009_ShardLongReads.unmapped_shards[main_shard_index]

        String fbmrq_prefix = basename(sharded_reads, ".bam")

        # Filter out the kinetics tags from PB files:
        call PB.RemoveKineticsTags as t_010_RemoveKineticsTags {
            input:
                bam = sharded_reads,
                prefix = SM + "_kinetics_removed"
        }

        # Handle setting up the things that we need for further processing of CCS-only reads:
        call PB.FindCCSReport as t_011_FindCCSReport {
            input:
                gcs_input_dir = gcs_input_dir
        }

        # 1 - filter the reads by the minimum read quality:
        call Utils.Bamtools as t_012_FilterS2EByMinReadQuality {
            input:
                bamfile = t_010_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_good_reads",
                cmd = "filter",
                args = '-tag "rq":">=' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        # 1.5 - Get the "rejected" reads:
        call Utils.Bamtools as t_013_GetS2ECcsRejectedReads {
            input:
                bamfile = t_010_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_rejected_reads",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        #################################################################################################################
        # 2 - Get reads we can reclaim:
        call Utils.Bamtools as t_014_ExtractS2ECcsReclaimableReads {
            input:
                bamfile = t_010_RemoveKineticsTags.bam_file,
                prefix = fbmrq_prefix + "_reads_for_ccs_reclamation",
                cmd = "filter",
                args = '-tag "rq":"<' + min_read_quality + '" -length "<=' + max_reclamation_length + '"',
                runtime_attr_override = disable_preemption_runtime_attrs
        }

        if ( ! t_008_CheckForAnnotatedReads.bam_has_annotations ) {
            # 3: Longbow annotate ccs reads
            call LONGBOW.Annotate as t_015_AnnotateS2ECCSReads {
                input:
                    reads = t_012_FilterS2EByMinReadQuality.bam_out,
                    model = mas_seq_model
            }
            # 4: Longbow annotate reclaimable reads
            call LONGBOW.Annotate as t_016_AnnotateS2EReclaimableReads {
                input:
                    reads = t_014_ExtractS2ECcsReclaimableReads.bam_out,
                    model = mas_seq_model
            }
        }

        File annotated_S2E_ccs_file = if t_008_CheckForAnnotatedReads.bam_has_annotations then t_012_FilterS2EByMinReadQuality.bam_out else select_first([t_015_AnnotateS2ECCSReads.annotated_bam])
        File annotated_S2E_reclaimable_file = if t_008_CheckForAnnotatedReads.bam_has_annotations then t_014_ExtractS2ECcsReclaimableReads.bam_out else select_first([t_016_AnnotateS2EReclaimableReads.annotated_bam])

        # 5: Longbow filter ccs annotated reads
        call LONGBOW.Filter as t_017_FilterS2ECCSReads {
            input:
                bam = annotated_S2E_ccs_file,
                prefix = SM + "_subshard",
                model = mas_seq_model
        }

        # 6: Longbow filter ccs reclaimable reads
        call LONGBOW.Filter as t_018_FilterS2EReclaimableReads {
            input:
                bam = annotated_S2E_reclaimable_file,
                prefix = SM + "_subshard",
                model = mas_seq_model
        }

        # 7: PBIndex CCS reads
        call PB.PBIndex as t_019_PbIndexS2ELongbowPassedCcsReads {
            input:
                bam = t_017_FilterS2ECCSReads.passed_reads
        }

        call PB.PBIndex as t_020_PbIndexS2ELongbowFailedCcsReads {
            input:
                bam = t_017_FilterS2ECCSReads.failed_reads
        }

        # 8: PBIndex reclaimable reads
        call PB.PBIndex as t_021_PbIndexS2ELongbowPassedReclaimedReads {
            input:
                bam = t_018_FilterS2EReclaimableReads.passed_reads
        }
        call PB.PBIndex as t_022_PbIndexS2ELongbowFailedReclaimableReads {
            input:
                bam = t_018_FilterS2EReclaimableReads.failed_reads
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


        # Segment CCS reads into array elements:
        call LONGBOW.Segment as t_023_SegmentS2ECcsReads {
            input:
                annotated_reads = t_017_FilterS2ECCSReads.passed_reads,
                prefix = SM + "_ccs_array_elements_shard_" + main_shard_index,
                extra_args = "-b",
                model = mas_seq_model
        }

        # Now call the new Longbow Sift to remove individual reads that are incomplete / truncated / malformed:
        call LONGBOW.Sift as t_024_LongbowSiftCCSArrayElements {
            input:
                segmented_input_reads = t_023_SegmentS2ECcsReads.segmented_bam,
                model = mas_seq_model,
                prefix = SM + "_ccs_array_elements_annotated_sifted_" + main_shard_index
        }

        # Now that we've annotated the reads, we can pad the UMIs by a couple of bases to aid in the deduping:
        call LONGBOW.Pad as t_025_LongbowPadCCSArrayElementUMIs {
            input:
                reads = t_024_LongbowSiftCCSArrayElements.sifted_bam,
                model = mas_seq_model,
                tag_to_expand = "ZU",
                padding = ccs_umi_padding,
                prefix = SM + "_ccs_array_elements_annotated_umi_padded_shard_" + main_shard_index
        }

        call LONGBOW.Pad as t_026_LongbowPadCCSArrayElementCBCs {
            input:
                reads = t_025_LongbowPadCCSArrayElementUMIs.padded_tag_bam,
                model = mas_seq_model,
                tag_to_expand = "CR",
                new_tag_dest = expanded_cbc_tag,
                padding = ccs_cbc_padding,
                prefix = SM + "_ccs_array_elements_annotated_cbc_padded_shard_" + main_shard_index
        }

        # Now we should correct our barcodes based on the whitelist:
        call LONGBOW.Correct as t_027_LongbowCorrectCCSCorrectedArrayElementCBCs {
            input:
                reads = t_026_LongbowPadCCSArrayElementCBCs.padded_tag_bam,
                barcode_allow_list = cell_barcode_whitelist,
                model = mas_seq_model,
                ccs_lev_dist_threshold = ccs_lev_dist,
                clr_lev_dist_threshold = clr_lev_dist,
                prefix = SM + "_ccs_array_elements_annotated_padded_cbc_corrected_shard_" + main_shard_index,
                raw_barcode_tag = expanded_cbc_tag,
                corrected_barcode_tag = "CB",
        }

        call TENX.AdjustUmiSequenceWithAdapterAlignment as t_028_AdjustCCSUMIs {
            input:
                bam = t_027_LongbowCorrectCCSCorrectedArrayElementCBCs.corrected_barcodes_bam,
                short_read_umis = short_read_umis_tsv,
                prefix = SM + "_ccs_array_elements_annotated_padded_cbc_corrected_UMI_adjusted_shard_" + main_shard_index,
        }

        call LONGBOW.Extract as t_029_LongbowExtractCcsArrayElements {
            input:
                bam = t_028_AdjustCCSUMIs.output_bam,
                prefix = SM + "_ccs_array_elements_cbc_umi_padded_extracted_shard_" + main_shard_index,
        }

        #####################
        # CCS Reclaimed / CLR:

        # Segment Reclaimed reads into array elements:
        call LONGBOW.Segment as t_030_SegmentS2ECcsReclaimedReads {
            input:
                annotated_reads = t_018_FilterS2EReclaimableReads.passed_reads,
                prefix = SM + "_ccs_reclaimed_array_elements_shard_" + main_shard_index,
                extra_args = "-b",
                model = mas_seq_model
        }

        # Now call the new Longbow Sift to remove individual reads that are incomplete / truncated / malformed:
        call LONGBOW.Sift as t_031_LongbowSiftCcsReclaimedArrayElements {
            input:
                segmented_input_reads = t_030_SegmentS2ECcsReclaimedReads.segmented_bam,
                model = mas_seq_model,
                prefix = SM + "_ccs_reclaimed_array_elements_annotated_sifted_" + main_shard_index
        }

        # Now that we've annotated the reads, we can pad the UMIs by a couple of bases to aid in the deduping:
        call LONGBOW.Pad as t_032_LongbowPadCcsReclaimedArrayElementUMIs {
            input:
                reads = t_031_LongbowSiftCcsReclaimedArrayElements.sifted_bam,
                model = mas_seq_model,
                tag_to_expand = "ZU",
                padding = ccs_umi_padding,
                prefix = SM + "_ccs_reclaimed_array_elements_annotated_umi_padded_shard_" + main_shard_index,
        }

        call LONGBOW.Pad as t_033_LongbowPadCcsReclaimedArrayElementCBCs {
            input:
                reads = t_032_LongbowPadCcsReclaimedArrayElementUMIs.padded_tag_bam,
                model = mas_seq_model,
                tag_to_expand = "CR",
                new_tag_dest = expanded_cbc_tag,
                padding = ccs_cbc_padding,
                prefix = SM + "_ccs_reclaimed_array_elements_annotated_cbc_padded_shard_" + main_shard_index,
        }

        # Now we should correct our barcodes based on the whitelist:
        call LONGBOW.Correct as t_034_LongbowCorrectCcsReclaimedArrayElementCBCs {
            input:
                reads = t_033_LongbowPadCcsReclaimedArrayElementCBCs.padded_tag_bam,
                barcode_allow_list = cell_barcode_whitelist,
                model = mas_seq_model,
                ccs_lev_dist_threshold = ccs_lev_dist,
                clr_lev_dist_threshold = clr_lev_dist,
                prefix = SM + "_ccs_reclaimed_array_elements_annotated_padded_cbc_corrected_shard_" + main_shard_index,
                raw_barcode_tag = expanded_cbc_tag,
                corrected_barcode_tag = "CB",
        }

        call TENX.AdjustUmiSequenceWithAdapterAlignment as t_035_AdjustCcsReclaimedUMIs {
            input:
                bam = t_034_LongbowCorrectCcsReclaimedArrayElementCBCs.corrected_barcodes_bam,
                short_read_umis = short_read_umis_tsv,
                prefix = SM + "_ccs_reclaimed_array_elements_annotated_padded_cbc_corrected_UMI_adjusted_shard_" + main_shard_index,
        }

        call LONGBOW.Extract as t_036_LongbowExtractCcsReclaimedArrayElements {
            input:
                bam = t_035_AdjustCcsReclaimedUMIs.output_bam,
                prefix = SM + "_ccs_reclaimed_array_elements_cbc_umi_padded_extracted_shard_" + main_shard_index,
        }

        ###############

        # Now align the array elements with their respective alignment presets.

        # Align CCS reads to the genome:
        call AR.Minimap2 as t_037_AlignCCSArrayElementsToGenome {
            input:
                reads      = [ t_029_LongbowExtractCcsArrayElements.extracted_bam ],
                ref_fasta  = ref_fasta,
                tags_to_preserve = tags_to_preserve,
                map_preset = "splice:hq",
                prefix = SM + "_ccs_array_elements_extracted_aligned_shard_"+ main_shard_index,
                runtime_attr_override = object { mem_gb: 32 }
        }

        call LONGBOW.TagFix as t_038_LongbowTagfixAlignedCcsArrayElements {
            input:
                bam = t_037_AlignCCSArrayElementsToGenome.aligned_bam,
                prefix = SM + "_ccs_array_elements_extracted_aligned_tagfixed_shard_" + main_shard_index,
        }

        # Align Reclaimed reads to the genome:
        call AR.Minimap2 as t_039_AlignReclaimedArrayElementsToGenome {
            input:
                reads      = [ t_036_LongbowExtractCcsReclaimedArrayElements.extracted_bam ],
                ref_fasta  = ref_fasta,
                tags_to_preserve = tags_to_preserve,
                map_preset = "splice",
                prefix = SM + "_ccs_reclaimed_array_elements_extracted_aligned_shard_" + main_shard_index,
                runtime_attr_override = object { mem_gb: 32 }
        }

        call LONGBOW.TagFix as t_040_LongbowTagfixAlignedCcsReclaimedArrayElements {
            input:
                bam = t_039_AlignReclaimedArrayElementsToGenome.aligned_bam,
                prefix = SM + "_ccs_reclaimed_array_elements_extracted_aligned_tagfixed_shard_" + main_shard_index,
        }
    }

    # Merge the arrays:
    call Utils.MergeBams as t_041_MergeCCSLongbowAnnotatedArrayReads {
        input:
            bams = annotated_S2E_ccs_file,
            prefix = SM + "_ccs_array_reads_longbow_annotated"
    }
    call PB.PBIndex as t_042_PbIndexMergedCCSLongbowAnnotatedArrayReads { input: bam = t_041_MergeCCSLongbowAnnotatedArrayReads.merged_bam }

    call Utils.MergeBams as t_043_MergeCCSReclaimableLongbowAnnotatedArrayReads {
        input:
            bams = annotated_S2E_reclaimable_file,
            prefix = SM + "_ccs_reclaimable_array_reads_longbow_annotated"
    }
    call PB.PBIndex as t_044_PbIndexMergedCCSReclaimableLongbowAnnotatedArrayReads { input: bam = t_043_MergeCCSReclaimableLongbowAnnotatedArrayReads.merged_bam }

    call Utils.MergeBams as t_045_MergeCCSLongbowPassedArrayReads {
        input:
            bams = t_017_FilterS2ECCSReads.passed_reads,
            prefix = SM + "_ccs_array_reads_longbow_passed"
    }
    call PB.PBIndex as t_046_PbIndexMergedCCSLongbowPassedArrayReads { input: bam = t_045_MergeCCSLongbowPassedArrayReads.merged_bam }

    call Utils.MergeBams as t_047_MergeCCSLongbowFailedArrayReads {
        input:
            bams = t_017_FilterS2ECCSReads.failed_reads,
            prefix = SM + "_ccs_array_reads_longbow_failed"
    }
    call PB.PBIndex as t_048_PbIndexMergedCCSLongbowFailedArrayReads { input: bam = t_047_MergeCCSLongbowFailedArrayReads.merged_bam }

    call Utils.MergeBams as t_049_MergeCCSReclaimedArrayReads {
        input:
            bams = t_018_FilterS2EReclaimableReads.passed_reads,
            prefix = SM + "_ccs_reclaimed_array_reads_longbow_passed"
    }
    call PB.PBIndex as t_050_PbIndexMergedCCSReclaimedArrayReads { input: bam = t_049_MergeCCSReclaimedArrayReads.merged_bam }

    call Utils.MergeBams as t_051_MergeCCSUnreclaimableArrayReads {
        input:
            bams = t_018_FilterS2EReclaimableReads.failed_reads,
            prefix = SM + "_ccs_unreclaimable_array_reads_longbow_failed"
    }
    call PB.PBIndex as t_052_PbIndexMergedCCSUnreclaimableReclaimedArrayReads { input: bam = t_051_MergeCCSUnreclaimableArrayReads.merged_bam }

    call Utils.MergeBams as t_053_MergeLongbowPassedReads {
        input:
            bams = flatten([t_017_FilterS2ECCSReads.passed_reads, t_018_FilterS2EReclaimableReads.passed_reads]),
            prefix = SM + "_longbow_passed_array_reads"
    }
    call PB.PBIndex as t_054_PbIndexMergedLongbowPassingReads { input: bam = t_053_MergeLongbowPassedReads.merged_bam }

    call Utils.MergeBams as t_055_MergeLongbowFailedReads {
        input:
            bams = flatten([t_017_FilterS2ECCSReads.failed_reads, t_018_FilterS2EReclaimableReads.failed_reads]),
            prefix = SM + "_longbow_failed_array_reads"
    }
    call PB.PBIndex as t_056_PbIndexMergedLongbowFailedReads { input: bam = t_055_MergeLongbowFailedReads.merged_bam }

    call Utils.MergeBams as t_057_MergeAllLongbowAnnotatedReads {
        input:
            bams = flatten([annotated_S2E_ccs_file, annotated_S2E_reclaimable_file]),
            prefix = SM + "_all_longbow_annotated_array_reads"
    }
    call PB.PBIndex as t_058_PbIndexMergedAllLongbowAnnotatedReads { input: bam = t_057_MergeAllLongbowAnnotatedReads.merged_bam }

    # Merge Filtered CCS reads together:
    call Utils.MergeBams as t_059_MergeCCSArrayElements {
        input:
            bams = t_023_SegmentS2ECcsReads.segmented_bam,
            prefix = SM + "_ccs_array_elements"
    }

    # Merge CCS Barcode Conf files:
    call Utils.MergeFiles as t_060_MergeCCSBarcodeConfShards {
        input:
            files_to_merge = t_023_SegmentS2ECcsReads.barcode_conf_file,
            merged_file_name = SM + "_ccs_array_element_barcode_confs.txt"
    }

    # Merge sifted CCS Array Elements:
    call Utils.MergeBams as t_061_MergeCCSArrayElementsSifted {
        input:
            bams = t_024_LongbowSiftCCSArrayElements.sifted_bam,
            prefix = SM + "_ccs_array_elements_annotated_sifted"
    }
    call Utils.MergeBams as t_062_MergeCCSArrayElementsSiftedFailed {
        input:
            bams = t_024_LongbowSiftCCSArrayElements.sift_failed_bam,
            prefix = SM + "_ccs_array_elements_annotated_sifted_failed"
    }

    # Merge UMI-Padded CCS Array Elements:
    call Utils.MergeBams as t_063_MergeCCSArrayElementsUmiPadded {
        input:
            bams = t_025_LongbowPadCCSArrayElementUMIs.padded_tag_bam,
            prefix = SM + "_ccs_array_elements_no_ends_umi_padded"
    }

    # Merge CBC-UMI-Padded CCS Array Elements:
    call Utils.MergeBams as t_064_MergeCCSArrayElementsUmiCbcPadded {
        input:
            bams = t_026_LongbowPadCCSArrayElementCBCs.padded_tag_bam,
            prefix = SM + "_ccs_array_elements_no_ends_cbc_umi_padded"
    }

    # Merge Corrected CBC CCS Array Elements:
    call Utils.MergeBams as t_065_MergeCCSArrayElementsUmiCbcPaddedCbcCorrected {
        input:
            bams = t_027_LongbowCorrectCCSCorrectedArrayElementCBCs.corrected_barcodes_bam,
            prefix = SM + "_ccs_array_elements_no_ends_cbc_umi_padded"
    }

    call Utils.MergeBams as t_066_MergeLongbowPaddedCBCUncorrectableCCSArrayElements {
        input:
            bams = t_027_LongbowCorrectCCSCorrectedArrayElementCBCs.uncorrected_barcodes_bam,
            prefix = SM + "_ccs_array_elements_aligned_annotated_padded_CBC_uncorrectable"
    }

    call Utils.MergeBams as t_067_MergeLongbowPaddedCBCCorrectedCCSArrayElements {
        input:
            bams = t_028_AdjustCCSUMIs.output_bam,
            prefix = SM + "_ccs_array_elements_aligned_annotated_padded_CBC_corrected"
    }

    call Utils.MergeBams as t_068_MergeLongbowExtractedCcsArrayElements {
        input:
            bams = t_029_LongbowExtractCcsArrayElements.extracted_bam,
            prefix = SM + "_ccs_array_elements_cbc_umi_padded_extracted"
    }

     # Merge Filtered CCS Reclaimed reads together:
    call Utils.MergeBams as t_069_MergeCCSReclaimedArrayElements {
        input:
            bams = t_030_SegmentS2ECcsReclaimedReads.segmented_bam,
            prefix = SM + "_ccs_reclaimed_array_elements"
    }

    # Merge CCS Barcode Conf files:
    call Utils.MergeFiles as t_070_MergeCCSReclaimedBarcodeConfShards {
        input:
            files_to_merge = t_030_SegmentS2ECcsReclaimedReads.barcode_conf_file,
            merged_file_name = SM + "_ccs_reclaimed_array_element_barcode_confs"
    }

    # Merge sifted CCS Reclaiomed Array Elements:
    call Utils.MergeBams as t_071_MergeCCSReclaimedArrayElementsSifted {
        input:
            bams = t_031_LongbowSiftCcsReclaimedArrayElements.sifted_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_annotated_sifted"
    }
    call Utils.MergeBams as t_072_MergeCCSReclaimedArrayElementsSiftedFailed {
        input:
            bams = t_031_LongbowSiftCcsReclaimedArrayElements.sift_failed_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_annotated_sifted_failed"
    }

    # Merge UMI-Padded CCS Array Elements:
    call Utils.MergeBams as t_073_MergeCCSReclaimedArrayElementsUmiPadded {
        input:
            bams = t_032_LongbowPadCcsReclaimedArrayElementUMIs.padded_tag_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_no_ends_umi_padded"
    }

    # Merge CBC-UMI-Padded CCS Array Elements:
    call Utils.MergeBams as t_074_MergeCCSReclaimedArrayElementsUmiCbcPadded {
        input:
            bams = t_033_LongbowPadCcsReclaimedArrayElementCBCs.padded_tag_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_no_ends_cbc_umi_padded"
    }

    # Merge Corrected CBC CCS Array Elements:
    call Utils.MergeBams as t_075_MergeCCSReclaimedArrayElementsUmiCbcPaddedCbcCorrected {
        input:
            bams = t_034_LongbowCorrectCcsReclaimedArrayElementCBCs.corrected_barcodes_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_no_ends_cbc_umi_padded"
    }

    call Utils.MergeBams as t_076_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElements {
        input:
            bams = t_034_LongbowCorrectCcsReclaimedArrayElementCBCs.uncorrected_barcodes_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_aligned_annotated_padded_CBC_uncorrectable"
    }

    call Utils.MergeBams as t_077_MergeLongbowPaddedCBCCorrectedCCSReclaimedArrayElementsShards {
        input:
            bams = t_035_AdjustCcsReclaimedUMIs.output_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_aligned_annotated_padded_CBC_corrected"
    }

    call Utils.MergeBams as t_078_MergeLongbowExtractedCcsReclaimedArrayElements {
        input:
            bams = t_036_LongbowExtractCcsReclaimedArrayElements.extracted_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_cbc_umi_padded_extracted"
    }

    call Utils.MergeBams as t_079_MergeAllRawArrayElements {
        input:
            bams = flatten([t_023_SegmentS2ECcsReads.segmented_bam, t_030_SegmentS2ECcsReclaimedReads.segmented_bam]),
            prefix = SM + "_raw_array_elements"
    }

    # Merge Aligned CCS array elements together:
    call Utils.MergeBams as t_080_MergeAlignedCCSArrayElements {
        input:
            bams = t_038_LongbowTagfixAlignedCcsArrayElements.tag_fixed_bam,
            prefix = SM + "_ccs_array_elements_padded_aligned_tagfixed"
    }

    # Merge Aligned CCS Reclaimed array elements together:
    call Utils.MergeBams as t_081_MergeAlignedCCSReclaimedArrayElements {
        input:
            bams = t_040_LongbowTagfixAlignedCcsReclaimedArrayElements.tag_fixed_bam,
            prefix = SM + "_ccs_reclaimed_array_elements_padded_aligned_tagfixed"
    }

    call Utils.MergeBams as t_082_MergeAllAlignedArrayElements {
        input:
            bams = flatten([t_038_LongbowTagfixAlignedCcsArrayElements.tag_fixed_bam, t_040_LongbowTagfixAlignedCcsReclaimedArrayElements.tag_fixed_bam]),
            prefix = SM + "_array_elements_padded_aligned_tagfixed"
    }

    # Merge CCS Barcode Conf files:
    call Utils.MergeFiles as t_083_MergeAllCCSBarcodeConfShards {
        input:
            files_to_merge = t_023_SegmentS2ECcsReads.segmented_bam,
            merged_file_name = SM + "_ccs_array_element_barcode_confs.txt"
    }

    # Merge CCS Barcode Conf files:
    call Utils.MergeFiles as t_084_MergeAllCCSReclaimedBarcodeConfShards {
        input:
            files_to_merge = t_030_SegmentS2ECcsReclaimedReads.segmented_bam,
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
    call Utils.FilterMasSeqReadsWithGatk as t_085_AlignmentFilterForCcsArrayElements {
        input:
            bam_file = t_080_MergeAlignedCCSArrayElements.merged_bam,
            bam_index = t_080_MergeAlignedCCSArrayElements.merged_bai,
            prefix = SM + "_CCS_ArrayElements_Annotated_Aligned_PrimaryOnly",
            runtime_attr_override = filterReadsAttrs
    }

    call Utils.FilterMasSeqReadsWithGatk as t_086_AlignmentFilterForReclaimedArrayElements {
        input:
            bam_file = t_081_MergeAlignedCCSReclaimedArrayElements.merged_bam,
            bam_index = t_081_MergeAlignedCCSReclaimedArrayElements.merged_bai,
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
    call Utils.MergeBams as t_087_MergeAllAlignedAndFilteredArrayElements {
        input:
            bams = [t_085_AlignmentFilterForCcsArrayElements.bam, t_086_AlignmentFilterForReclaimedArrayElements.bam],
            prefix = SM + "_all_array_elements_aligned_for_txome_discovery"
    }

    call StringTie2.Quantify as t_088_ST2_Quant {
        input:
            aligned_bam = t_087_MergeAllAlignedAndFilteredArrayElements.merged_bam,
            aligned_bai = t_087_MergeAllAlignedAndFilteredArrayElements.merged_bai,
            gtf = genome_annotation_gtf,
            keep_retained_introns = false,
            prefix = SM + "_StringTie2_Quantify",
    }

    call StringTie2.ExtractTranscriptSequences as t_089_ST2_ExtractTranscriptSequences  {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_index,
            gtf = t_088_ST2_Quant.st_gtf,
            prefix = SM + "_StringTie2_ExtractTranscriptSequences",
    }

    call StringTie2.CompareTranscriptomes as t_090_ST2_CompareTranscriptomes {
        input:
            guide_gtf = genome_annotation_gtf,
            new_gtf = t_088_ST2_Quant.st_gtf,
            prefix = SM + "_StringTie2_CompareTranscriptome",
    }

    ##########################################################################################################################
    ##########################################################################################################################

    #################
    # Here we restore the original read names to the bam because we're hashing them with Longbow.segment:

    # Restore original read names to CCS reads:
    call TX_PRE.RestoreOriginalReadNames as t_091_RestoreCcsOriginalReadNames {
        input:
            bam = t_085_AlignmentFilterForCcsArrayElements.bam,
            prefix =  SM + "_CCS_cbc_annotated_array_elements_padded_original_names"
    }

    # Restore original read names to CLR reads:
    call TX_PRE.RestoreOriginalReadNames as t_092_RestoreClrOriginalReadNames {
        input:
            bam = t_086_AlignmentFilterForReclaimedArrayElements.bam,
            prefix =  SM + "_CLR_cbc_annotated_array_elements_padded_original_names"
    }

    # Merge Aligned CCS and Reclaimed reads together:
    call Utils.MergeBams as t_093_MergeAllAnnotatedArrayElementsWithOriginalNames {
        input:
            bams = [t_091_RestoreCcsOriginalReadNames.bam_out, t_092_RestoreClrOriginalReadNames.bam_out],
            prefix = SM + "_all_cbc_annotated_array_elements_padded_original_names"
    }

    #################
    # Now we have to split the reads again, process them into gff files, run gffcompare and then aggregate the results in a graph

    # We can actually compare the references without needing to scatter:
    call TX_PRE.GffCompare as t_094_GffCompareStringtie2toGencode {
        input:
            gff_ref = t_088_ST2_Quant.st_gtf,
            gff_query = genome_annotation_gtf,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
    }
    call TX_PRE.GffCompare as t_095_GffCompareGencodetoStringtie2 {
        input:
            gff_ref = genome_annotation_gtf,
            gff_query = t_088_ST2_Quant.st_gtf,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
    }

    # Split by contig:
    call TX_PRE.SplitBamByContig as t_096_SplitArrayElementsByContig {
        input:
            bam = t_093_MergeAllAnnotatedArrayElementsWithOriginalNames.merged_bam,
            prefix = SM + "_all_cbc_annotated_array_elements_padded_original_names"
    }

    # For each contig:
    scatter (i in range(length(t_096_SplitArrayElementsByContig.contig_bams))) {

        File contig_bam = t_096_SplitArrayElementsByContig.contig_bams[i]
        String contig_name = t_096_SplitArrayElementsByContig.contig_names[i]

        # Create a GFF file:
        call TX_PRE.ConvertSplicedBamToGff as t_097_ConvertSplicedBamToGff {
            input:
                bam = contig_bam
        }

        # Compare GFF files:
        call TX_PRE.GffCompare as t_098_GffCompareStringtie2toMasSeqReads {
            input:
                gff_ref = t_088_ST2_Quant.st_gtf,
                gff_query = t_097_ConvertSplicedBamToGff.gff,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
        }

        call TX_PRE.GffCompare as t_099_GffCompareGencodetoMasSeqReads {
            input:
                gff_ref = genome_annotation_gtf,
                gff_query = t_097_ConvertSplicedBamToGff.gff,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
        }

        # Create the comparison graph and tsv files:
        call TX_POST.QuantifyGffComparison as t_100_QuantifyGffComparison {
            input:
                genome_gtf = genome_annotation_gtf,
                st2_gencode_refmap = t_094_GffCompareStringtie2toGencode.refmap,
                st2_gencode_tmap = t_094_GffCompareStringtie2toGencode.tmap,
                st2_read_refmap = t_098_GffCompareStringtie2toMasSeqReads.refmap,
                st2_read_tmap = t_098_GffCompareStringtie2toMasSeqReads.tmap,
                gencode_st2_refmap = t_095_GffCompareGencodetoStringtie2.refmap,
                gencode_st2_tmap = t_095_GffCompareGencodetoStringtie2.tmap,
                gencode_read_refmap = t_099_GffCompareGencodetoMasSeqReads.refmap,
                gencode_read_tmap = t_099_GffCompareGencodetoMasSeqReads.tmap,
                prefix = SM + "_all_cbc_annotated_array_elements_padded_" + contig_name
        }
    }

    # Merge our tx equivalance classes assignments and eq classes:
    call TX_POST.CombineEqClassFiles as t_101_CombineEqClassFiles {
        input:
            gene_eq_class_definitions = t_100_QuantifyGffComparison.gene_eq_class_labels_file,
            gene_assignment_files = t_100_QuantifyGffComparison.gene_assignments_file,
            equivalence_class_definitions = t_100_QuantifyGffComparison.tx_equivalence_class_labels_file,
            equivalence_classes = t_100_QuantifyGffComparison.tx_equivalence_class_file,
            prefix = SM + "_all_cbc_annotated_array_elements_padded"
    }

    ############################################################
    # Quantify Transcripts:
    ##########

    # Use old quant method here as a baseline for comparison:
    call TX_POST.CopyEqClassInfoToTag as t_102_CopyEqClassInfoToTag {
        input:
            bam = t_093_MergeAllAnnotatedArrayElementsWithOriginalNames.merged_bam,
            eq_class_file = t_101_CombineEqClassFiles.combined_tx_eq_class_assignments,
            prefix = SM + "_annotated_array_elements_for_quant_with_gene_names"
    }

    call LONGBOW.Correct_UMI as t_103_LongbowCorrectUmi {
        input:
            bam = t_102_CopyEqClassInfoToTag.bam_out,
            prefix = SM + "_annotated_array_elements_for_quant_with_gene_names"
    }

    # Because of how we're doing things, we need to pull out the CCS and CCS Reclaimed reads from the output of the
    # set cover correction:
    call Utils.Bamtools as t_104_GetCcsCorrectedReadsWithCorrectedUmis {
        input:
            bamfile = t_103_LongbowCorrectUmi.umi_corrected_bam,
            prefix = SM + "_annotated_array_elements_for_quant_with_gene_names.corrected_umis.CCS",
            cmd = "filter",
            args = '-tag "rq":">=' + min_read_quality + '"',
            runtime_attr_override = disable_preemption_runtime_attrs
    }
    call Utils.IndexBam as t_105_IndexCcsReadsWithCorrectedUmis {input: bam = t_104_GetCcsCorrectedReadsWithCorrectedUmis.bam_out }

    call Utils.Bamtools as t_106_GetCcsReclaimedReadsWithCorrectedUmis {
        input:
            bamfile = t_103_LongbowCorrectUmi.umi_corrected_bam,
            prefix = SM + "_annotated_array_elements_for_quant_with_gene_names.corrected_umis.CCS_Reclaimed",
            cmd = "filter",
            args = '-tag "rq":"<' + min_read_quality + '"',
            runtime_attr_override = disable_preemption_runtime_attrs
    }
    call Utils.IndexBam as t_107_IndexCcsReclaimedReadsWithCorrectedUmis {input: bam = t_106_GetCcsReclaimedReadsWithCorrectedUmis.bam_out }

    call UMI_TOOLS.Run_Group as t_108_UMIToolsGroup {
        input:
            aligned_transcriptome_reads = t_103_LongbowCorrectUmi.umi_corrected_bam,
            aligned_transcriptome_reads_index = t_103_LongbowCorrectUmi.umi_corrected_bam_index,
            do_per_cell = true,
            prefix = SM + "_annotated_array_elements_with_gene_names_with_umi_tools_group_correction"
    }

    # Create CCS count matrix and anndata:
    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_109_CreateCCSCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_104_GetCcsCorrectedReadsWithCorrectedUmis.bam_out,
            tx_equivalence_class_assignments = t_101_CombineEqClassFiles.combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = SM + "_ccs_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_110_CreateCCSCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_109_CreateCCSCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = t_088_ST2_Quant.st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = t_101_CombineEqClassFiles.combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = t_101_CombineEqClassFiles.combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = t_101_CombineEqClassFiles.combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = t_101_CombineEqClassFiles.combined_gene_eq_class_defs,
            prefix = SM + "_ccs_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: 64}
    }

    # Create CLR count matrix and anndata:
    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_111_CreateCLRCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_106_GetCcsReclaimedReadsWithCorrectedUmis.bam_out,
            tx_equivalence_class_assignments = t_101_CombineEqClassFiles.combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = SM + "_clr_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_112_CreateCLRCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_111_CreateCLRCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = t_088_ST2_Quant.st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = t_101_CombineEqClassFiles.combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = t_101_CombineEqClassFiles.combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = t_101_CombineEqClassFiles.combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = t_101_CombineEqClassFiles.combined_gene_eq_class_defs,
            prefix = SM + "_clr_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: 64}
    }

    # Create overall count matrix and anndata:
    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_113_CreateOverallCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_103_LongbowCorrectUmi.umi_corrected_bam,
            tx_equivalence_class_assignments = t_101_CombineEqClassFiles.combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = SM + "_overall_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_114_CreateOverallCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_113_CreateOverallCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = t_088_ST2_Quant.st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = t_101_CombineEqClassFiles.combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = t_101_CombineEqClassFiles.combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = t_101_CombineEqClassFiles.combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = t_101_CombineEqClassFiles.combined_gene_eq_class_defs,
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

    call AM.SamtoolsStats as t_115_BaselineArrayElementStats {
        input:
            bam = t_079_MergeAllRawArrayElements.merged_bam
    }

    call AM.SamtoolsStats as t_116_AlignedArrayElementStats {
        input:
            bam = t_082_MergeAllAlignedArrayElements.merged_bam
    }

    call AM.SamtoolsStats as t_117_AlignedFilteredArrayElementStats {
        input:
            bam = t_087_MergeAllAlignedAndFilteredArrayElements.merged_bam
    }

    call AM.SamtoolsStats as t_118_AlignedAnnotatedArrayElementsForQuantStats {
        input:
            bam = t_103_LongbowCorrectUmi.umi_corrected_bam
    }

    call LONGBOW.AggregateCorrectLogStats as t_119_AggregateLongbowCorrectStats {
        input:
            longbow_correct_log_files = flatten([t_034_LongbowCorrectCcsReclaimedArrayElementCBCs.log, t_027_LongbowCorrectCCSCorrectedArrayElementCBCs.log]),
            out_name = SM + "_longbow_correct_stats.txt"
    }

    # Get stats on CCS reads:
    call LONGBOW.Stats as t_120_CCS_longbow_stats {
        input:
            reads = t_041_MergeCCSLongbowAnnotatedArrayReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_CCS_Corrected",
    }

    # Get stats on Reclaimable reads:
    call LONGBOW.Stats as t_121_Reclaimable_longbow_stats {
        input:
            reads = t_043_MergeCCSReclaimableLongbowAnnotatedArrayReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_CCS_Reclaimable",
    }

    # Get stats on Reclaimed reads:
    call LONGBOW.Stats as t_122_Reclaimed_longbow_stats {
        input:
            reads = t_049_MergeCCSReclaimedArrayReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_CCS_Reclaimed",
    }

    # Get stats on All Passing reads (overall stats):
    call LONGBOW.Stats as t_123_Passed_longbow_stats {
        input:
            reads = t_053_MergeLongbowPassedReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_All_Longbow_Passed",
    }

    # Get stats on All Failed reads (overall stats):
    call LONGBOW.Stats as t_124_Failed_longbow_stats {
        input:
            reads = t_055_MergeLongbowFailedReads.merged_bam,
            model = mas_seq_model,
            prefix = SM + "_All_Longbow_Failed",
    }

    # Get stats on All reads (overall stats):
    call LONGBOW.Stats as t_125_Overall_longbow_stats {
        input:
            reads = t_057_MergeAllLongbowAnnotatedReads.merged_bam,
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

#    File keyfile = t_114_CreateOverallCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad

    # This seems to take longer to get to:
    File keyfile = t_108_UMIToolsGroup.output_tsv

    String base_out_dir = outdir + "/" + DIR + out_dir_suffix + "/" + t_001_WdlExecutionStartTimestamp.timestamp_string
    String stats_out_dir = base_out_dir + "/stats"
    String array_element_dir = base_out_dir + "/annotated_array_elements"
    String intermediate_reads_dir = base_out_dir + "/intermediate_reads"

    String meta_files_dir = base_out_dir + "/meta_files"

    String intermediate_array_reads_dir = intermediate_reads_dir + "/array_reads"
    String intermediate_array_elements_dir = intermediate_reads_dir + "/array_elements"

    String quant_dir = base_out_dir + "/quant"

    ##############################################################################################################
    # Finalize gene / tx assignments:
    call FF.FinalizeToDir as t_126_FinalizeEqClasses {
        input:
            files = [
                t_101_CombineEqClassFiles.combined_gene_eq_class_defs,
                t_101_CombineEqClassFiles.combined_gene_eq_class_assignments,
                t_101_CombineEqClassFiles.combined_tx_eq_class_defs,
                t_101_CombineEqClassFiles.combined_tx_eq_class_assignments,
            ],
            outdir = quant_dir + "/eqivalence_classes",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_127_FinalizeUmiToolsOutputs {
        input:
            files = [
                t_108_UMIToolsGroup.output_bam,
                t_108_UMIToolsGroup.output_tsv,
            ],
            outdir = quant_dir + "/UMITools",
            keyfile = keyfile
    }

    # CCS:
    call FF.FinalizeToDir as t_128_FinalizeCCSTxAndGeneAssignments {
        input:
            files = [
                t_109_CreateCCSCountMatrixFromAnnotatedBam.count_matrix,
                t_110_CreateCCSCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad,
            ],
            outdir = quant_dir + "/CCS",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_129_FinalizeCCSRawQuantPickles {
        input:
            files = t_110_CreateCCSCountMatrixAnndataFromEqClasses.pickles,
            outdir = quant_dir + "/CCS",
            keyfile = keyfile
    }

    # CLR:
    call FF.FinalizeToDir as t_130_FinalizeCLRTxAndGeneAssignments {
        input:
            files = [
                t_111_CreateCLRCountMatrixFromAnnotatedBam.count_matrix,
                t_112_CreateCLRCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad,
            ],
            outdir = quant_dir + "/CLR",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_131_FinalizeCLRRawQuantPickles {
        input:
            files = t_112_CreateCLRCountMatrixAnndataFromEqClasses.pickles,
            outdir = quant_dir + "/CLR",
            keyfile = keyfile
    }

    # Overall:
    call FF.FinalizeToDir as t_132_FinalizeOverallTxAndGeneAssignments {
        input:
            files = [
                t_113_CreateOverallCountMatrixFromAnnotatedBam.count_matrix,
                t_114_CreateOverallCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad,
            ],
            outdir = quant_dir + "/Overall",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_133_FinalizeOverallRawQuantPickles {
        input:
            files = t_114_CreateOverallCountMatrixAnndataFromEqClasses.pickles,
            outdir = quant_dir + "/Overall",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_134_FinalizeRefAndSt2Comparisons {
        input:
            files = [
                t_094_GffCompareStringtie2toGencode.refmap,
                t_094_GffCompareStringtie2toGencode.tmap,
                t_094_GffCompareStringtie2toGencode.tracking,
                t_094_GffCompareStringtie2toGencode.loci,
                t_094_GffCompareStringtie2toGencode.annotated_gtf,
                t_094_GffCompareStringtie2toGencode.stats,
                t_094_GffCompareStringtie2toGencode.log,

                t_095_GffCompareGencodetoStringtie2.refmap,
                t_095_GffCompareGencodetoStringtie2.tmap,
                t_095_GffCompareGencodetoStringtie2.tracking,
                t_095_GffCompareGencodetoStringtie2.loci,
                t_095_GffCompareGencodetoStringtie2.annotated_gtf,
                t_095_GffCompareGencodetoStringtie2.stats,
                t_095_GffCompareGencodetoStringtie2.log,
            ],
            outdir = quant_dir + "/gencode_and_stringtie2",
            keyfile = keyfile
    }

    # Finalize gene / tx assignment by contig:
    # NOTE: According to the scatter/gather documentation in the WDL spec, this will work correctly
    #       (https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#scatter--gather)
    scatter (i in range(length(t_096_SplitArrayElementsByContig.contig_bams))) {
        String contig = t_096_SplitArrayElementsByContig.contig_names[i]

        call FF.FinalizeToDir as t_135_FinalizeTxAndGeneAssignmentsByContig {
            input:
                files = [
                    t_098_GffCompareStringtie2toMasSeqReads.refmap[i],
                    t_098_GffCompareStringtie2toMasSeqReads.tmap[i],
                    t_098_GffCompareStringtie2toMasSeqReads.tracking[i],
                    t_098_GffCompareStringtie2toMasSeqReads.loci[i],
                    t_098_GffCompareStringtie2toMasSeqReads.annotated_gtf[i],
                    t_098_GffCompareStringtie2toMasSeqReads.stats[i],
                    t_098_GffCompareStringtie2toMasSeqReads.log[i],

                    t_099_GffCompareGencodetoMasSeqReads.refmap[i],
                    t_099_GffCompareGencodetoMasSeqReads.tmap[i],
                    t_099_GffCompareGencodetoMasSeqReads.tracking[i],
                    t_099_GffCompareGencodetoMasSeqReads.loci[i],
                    t_099_GffCompareGencodetoMasSeqReads.annotated_gtf[i],
                    t_099_GffCompareGencodetoMasSeqReads.stats[i],
                    t_099_GffCompareGencodetoMasSeqReads.log[i],

                    t_100_QuantifyGffComparison.gene_assignments_file[i],
                    t_100_QuantifyGffComparison.gene_eq_class_labels_file[i],
                    t_100_QuantifyGffComparison.tx_equivalence_class_labels_file[i],
                    t_100_QuantifyGffComparison.tx_equivalence_class_file[i],
                    t_100_QuantifyGffComparison.graph_gpickle[i],
                ],
                outdir = quant_dir + "/by_contig/" + contig,
                keyfile = keyfile
        }
    }

    ##############################################################################################################
    # Finalize annotated, aligned array elements:
    call FF.FinalizeToDir as t_136_FinalizeIntermediateAnnotatedArrayElements {
        input:
            files = [
                t_059_MergeCCSArrayElements.merged_bam,
                t_059_MergeCCSArrayElements.merged_bai,
                t_061_MergeCCSArrayElementsSifted.merged_bam,
                t_061_MergeCCSArrayElementsSifted.merged_bai,
                t_062_MergeCCSArrayElementsSiftedFailed.merged_bam,
                t_062_MergeCCSArrayElementsSiftedFailed.merged_bai,
                t_063_MergeCCSArrayElementsUmiPadded.merged_bam,
                t_063_MergeCCSArrayElementsUmiPadded.merged_bai,
                t_064_MergeCCSArrayElementsUmiCbcPadded.merged_bam,
                t_064_MergeCCSArrayElementsUmiCbcPadded.merged_bai,
                t_065_MergeCCSArrayElementsUmiCbcPaddedCbcCorrected.merged_bam,
                t_065_MergeCCSArrayElementsUmiCbcPaddedCbcCorrected.merged_bai,
                t_067_MergeLongbowPaddedCBCCorrectedCCSArrayElements.merged_bam,
                t_067_MergeLongbowPaddedCBCCorrectedCCSArrayElements.merged_bai,
                t_066_MergeLongbowPaddedCBCUncorrectableCCSArrayElements.merged_bam,
                t_066_MergeLongbowPaddedCBCUncorrectableCCSArrayElements.merged_bai,

                t_069_MergeCCSReclaimedArrayElements.merged_bam,
                t_069_MergeCCSReclaimedArrayElements.merged_bai,
                t_071_MergeCCSReclaimedArrayElementsSifted.merged_bam,
                t_071_MergeCCSReclaimedArrayElementsSifted.merged_bai,
                t_072_MergeCCSReclaimedArrayElementsSiftedFailed.merged_bam,
                t_072_MergeCCSReclaimedArrayElementsSiftedFailed.merged_bai,
                t_073_MergeCCSReclaimedArrayElementsUmiPadded.merged_bam,
                t_073_MergeCCSReclaimedArrayElementsUmiPadded.merged_bai,
                t_074_MergeCCSReclaimedArrayElementsUmiCbcPadded.merged_bam,
                t_074_MergeCCSReclaimedArrayElementsUmiCbcPadded.merged_bai,
                t_075_MergeCCSReclaimedArrayElementsUmiCbcPaddedCbcCorrected.merged_bam,
                t_075_MergeCCSReclaimedArrayElementsUmiCbcPaddedCbcCorrected.merged_bai,
                t_076_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElements.merged_bam,
                t_076_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElements.merged_bai,

                t_079_MergeAllRawArrayElements.merged_bam,
                t_079_MergeAllRawArrayElements.merged_bai,

                t_082_MergeAllAlignedArrayElements.merged_bam,
                t_082_MergeAllAlignedArrayElements.merged_bai,

                t_087_MergeAllAlignedAndFilteredArrayElements.merged_bam,
                t_087_MergeAllAlignedAndFilteredArrayElements.merged_bai,

                t_068_MergeLongbowExtractedCcsArrayElements.merged_bam,
                t_068_MergeLongbowExtractedCcsArrayElements.merged_bai,
                t_075_MergeCCSReclaimedArrayElementsUmiCbcPaddedCbcCorrected.merged_bam,
                t_075_MergeCCSReclaimedArrayElementsUmiCbcPaddedCbcCorrected.merged_bai,
                t_076_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElements.merged_bam,
                t_076_MergeLongbowPaddedCBCUncorrectableCCSReclaimedArrayElements.merged_bai,
                t_078_MergeLongbowExtractedCcsReclaimedArrayElements.merged_bam,
                t_078_MergeLongbowExtractedCcsReclaimedArrayElements.merged_bai,

                t_091_RestoreCcsOriginalReadNames.bam_out,
                t_092_RestoreClrOriginalReadNames.bam_out,

                t_093_MergeAllAnnotatedArrayElementsWithOriginalNames.merged_bam,
                t_093_MergeAllAnnotatedArrayElementsWithOriginalNames.merged_bai,
            ],
            outdir = intermediate_array_elements_dir,
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_137_FinalizeAnnotatedArrayElements {
        input:
            files = [
                t_103_LongbowCorrectUmi.umi_corrected_bam,
                t_103_LongbowCorrectUmi.umi_corrected_bam_index,

                t_104_GetCcsCorrectedReadsWithCorrectedUmis.bam_out,
                t_106_GetCcsReclaimedReadsWithCorrectedUmis.bam_out,

                t_087_MergeAllAlignedAndFilteredArrayElements.merged_bam,
                t_087_MergeAllAlignedAndFilteredArrayElements.merged_bai
            ],
            outdir = array_element_dir,
            keyfile = keyfile
    }

    call FF.FinalizeToFile as t_138_FinalizeCcsArrayElementCorrectedUmiIndex {
        input:
            file = t_105_IndexCcsReadsWithCorrectedUmis.bai,
            outfile = array_element_dir + "/" + SM + "_annotated_array_elements_for_quant_with_gene_names.corrected_umis.CCS.bam.bai",
            keyfile = keyfile
    }

    call FF.FinalizeToFile as t_139_FinalizeCcsReclaimedArrayElementCorrectedUmiIndex {
        input:
            file = t_107_IndexCcsReclaimedReadsWithCorrectedUmis.bai,
            outfile = array_element_dir + "/" + SM + "_annotated_array_elements_for_quant_with_gene_names.corrected_umis.CCS_Reclaimed",
            keyfile = keyfile
    }

    ##############################################################################################################
    # Finalize meta files:
    call FF.FinalizeToDir as t_140_FinalizeMeta {
        input:
            files = [
                cell_barcode_whitelist,
                t_083_MergeAllCCSBarcodeConfShards.merged_file,
                t_084_MergeAllCCSReclaimedBarcodeConfShards.merged_file
            ],
            outdir = meta_files_dir,
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_141_FinalizeCCSCBCcorrectionLogsToMeta {
        input:
            files = t_027_LongbowCorrectCCSCorrectedArrayElementCBCs.log,
            outdir = meta_files_dir + "/" + "ccs_cbc_correction_logs",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_142_FinalizeCCSRejectedCBCcorrectionLogsToMeta {
        input:
            files = t_034_LongbowCorrectCcsReclaimedArrayElementCBCs.log,
            outdir = meta_files_dir + "/" + "ccs_rejected_cbc_correction_logs",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_143_FinalizeCCSUmiAdjustmentLogs {
        input:
            files = t_028_AdjustCCSUMIs.log,
            outdir = meta_files_dir + "/" + "umi_adjustment_logs_ccs",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_144_FinalizeCCSReclaimedUmiAdjustmentLogs {
        input:
            files = t_035_AdjustCcsReclaimedUMIs.log,
            outdir = meta_files_dir + "/" + "umi_adjustment_logs_ccs_reclaimed",
            keyfile = keyfile
    }

    ##############################################################################################################
    # Finalize the discovered transcriptome:
    if ( !is_SIRV_data ) {
        call FF.FinalizeToDir as t_145_FinalizeDiscoveredTranscriptome {
            input:
                files = [
                    t_088_ST2_Quant.st_gtf,
                    t_089_ST2_ExtractTranscriptSequences.transcripts_fa,
                    t_089_ST2_ExtractTranscriptSequences.transcripts_fai,
                    t_089_ST2_ExtractTranscriptSequences.transcripts_dict,
                    t_090_ST2_CompareTranscriptomes.annotated_gtf,
                    t_090_ST2_CompareTranscriptomes.loci,
                    t_090_ST2_CompareTranscriptomes.stats,
                    t_090_ST2_CompareTranscriptomes.tracking,
                    t_090_ST2_CompareTranscriptomes.refmap,
                    t_090_ST2_CompareTranscriptomes.tmap,
                ],
                outdir = base_out_dir + "/discovered_transcriptome",
                keyfile = keyfile
        }
    }
    ##############################################################################################################
    # Finalize the intermediate reads files (from raw CCS corrected reads through split array elements)
    call FF.FinalizeToDir as t_146_FinalizeArrayReads {
        input:
            files = [

                t_041_MergeCCSLongbowAnnotatedArrayReads.merged_bam,
                t_041_MergeCCSLongbowAnnotatedArrayReads.merged_bai,
                t_042_PbIndexMergedCCSLongbowAnnotatedArrayReads.pbindex,

                t_043_MergeCCSReclaimableLongbowAnnotatedArrayReads.merged_bam,
                t_043_MergeCCSReclaimableLongbowAnnotatedArrayReads.merged_bai,
                t_044_PbIndexMergedCCSReclaimableLongbowAnnotatedArrayReads.pbindex,

                t_053_MergeLongbowPassedReads.merged_bam,
                t_053_MergeLongbowPassedReads.merged_bai,
                t_054_PbIndexMergedLongbowPassingReads.pbindex,

                t_055_MergeLongbowFailedReads.merged_bam,
                t_055_MergeLongbowFailedReads.merged_bai,
                t_056_PbIndexMergedLongbowFailedReads.pbindex,

                t_057_MergeAllLongbowAnnotatedReads.merged_bam,
                t_057_MergeAllLongbowAnnotatedReads.merged_bai,
                t_058_PbIndexMergedAllLongbowAnnotatedReads.pbindex,

                t_045_MergeCCSLongbowPassedArrayReads.merged_bam,
                t_045_MergeCCSLongbowPassedArrayReads.merged_bai,
                t_046_PbIndexMergedCCSLongbowPassedArrayReads.pbindex,

                t_047_MergeCCSLongbowFailedArrayReads.merged_bam,
                t_047_MergeCCSLongbowFailedArrayReads.merged_bai,
                t_048_PbIndexMergedCCSLongbowFailedArrayReads.pbindex,

                t_049_MergeCCSReclaimedArrayReads.merged_bam,
                t_049_MergeCCSReclaimedArrayReads.merged_bai,
                t_050_PbIndexMergedCCSReclaimedArrayReads.pbindex,

                t_051_MergeCCSUnreclaimableArrayReads.merged_bam,
                t_051_MergeCCSUnreclaimableArrayReads.merged_bai,
                t_052_PbIndexMergedCCSUnreclaimableReclaimedArrayReads.pbindex,

            ],
            outdir = intermediate_array_reads_dir,
            keyfile = keyfile
    }

    ##############################################################################################################
    # Finalize Stats:

    # Write out completion file so in the future we can be 100% sure that this run was good:
    call FF.FinalizeToDir as t_147_FinalizeHighLevelStats {
        input:
            files = [ t_011_FindCCSReport.ccs_report[0], t_119_AggregateLongbowCorrectStats.stats ],
            outdir = stats_out_dir,
            keyfile = keyfile
    }

    # Finalize Longbow Sift stats:
    scatter (i_f2 in range(length(t_009_ShardLongReads.unmapped_shards))) {
        call FF.FinalizeToFile as t_148_FinalizeCcsLongbowSiftStats {
            input:
                file = t_024_LongbowSiftCCSArrayElements.stats_tsv[i_f2],
                outfile = stats_out_dir + "/longbow_stats/sift/ccs/" + SM + "_ccs_array_elements_annotated_sifted_" + i_f2 + ".stats.tsv",
                keyfile = keyfile
        }
        call FF.FinalizeToFile as t_149_FinalizeCcsLongbowSiftSummaryStats {
            input:
                file = t_024_LongbowSiftCCSArrayElements.summary_stats_tsv[i_f2],
                outfile = stats_out_dir + "/longbow_stats/sift/ccs/" + SM + "_ccs_array_elements_annotated_sifted_" + i_f2 + ".summary_stats.tsv",
                keyfile = keyfile
        }
        call FF.FinalizeToFile as t_150_FinalizeCcsReclaimedLongbowSiftStats {
            input:
                file = t_031_LongbowSiftCcsReclaimedArrayElements.stats_tsv[i_f2],
                outfile = stats_out_dir + "/longbow_stats/sift/ccs_reclaimed/" + SM + "_ccs_reclaimed_array_elements_annotated_sifted_" + i_f2 + ".stats.tsv",
                keyfile = keyfile
        }
        call FF.FinalizeToFile as t_151_FinalizeCcsReclaimedLongbowSiftSummaryStats {
            input:
                file = t_031_LongbowSiftCcsReclaimedArrayElements.summary_stats_tsv[i_f2],
                outfile = stats_out_dir + "/longbow_stats/sift/ccs_reclaimed/" + SM + "_ccs_reclaimed_array_elements_annotated_sifted_" + i_f2 + ".summary_stats.tsv",
                keyfile = keyfile
        }
    }

    call FF.FinalizeToDir as t_152_FinalizeQuantArrayElementStats {
        input:
            files = [
                t_118_AlignedAnnotatedArrayElementsForQuantStats.raw_stats,
                t_118_AlignedAnnotatedArrayElementsForQuantStats.summary_stats,
                t_118_AlignedAnnotatedArrayElementsForQuantStats.first_frag_qual,
                t_118_AlignedAnnotatedArrayElementsForQuantStats.last_frag_qual,
                t_118_AlignedAnnotatedArrayElementsForQuantStats.first_frag_gc_content,
                t_118_AlignedAnnotatedArrayElementsForQuantStats.last_frag_gc_content,
                t_118_AlignedAnnotatedArrayElementsForQuantStats.acgt_content_per_cycle,
                t_118_AlignedAnnotatedArrayElementsForQuantStats.insert_size,
                t_118_AlignedAnnotatedArrayElementsForQuantStats.read_length_dist,
                t_118_AlignedAnnotatedArrayElementsForQuantStats.indel_distribution,
                t_118_AlignedAnnotatedArrayElementsForQuantStats.indels_per_cycle,
                t_118_AlignedAnnotatedArrayElementsForQuantStats.coverage_distribution,
                t_118_AlignedAnnotatedArrayElementsForQuantStats.gc_depth,
            ],
            outdir = stats_out_dir + "/array_elements_for_quant/",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_153_FinalizeTxomeDiscoveryArrayElementStats {
        input:
            files = [
                t_117_AlignedFilteredArrayElementStats.raw_stats,
                t_117_AlignedFilteredArrayElementStats.summary_stats,
                t_117_AlignedFilteredArrayElementStats.first_frag_qual,
                t_117_AlignedFilteredArrayElementStats.last_frag_qual,
                t_117_AlignedFilteredArrayElementStats.first_frag_gc_content,
                t_117_AlignedFilteredArrayElementStats.last_frag_gc_content,
                t_117_AlignedFilteredArrayElementStats.acgt_content_per_cycle,
                t_117_AlignedFilteredArrayElementStats.insert_size,
                t_117_AlignedFilteredArrayElementStats.read_length_dist,
                t_117_AlignedFilteredArrayElementStats.indel_distribution,
                t_117_AlignedFilteredArrayElementStats.indels_per_cycle,
                t_117_AlignedFilteredArrayElementStats.coverage_distribution,
                t_117_AlignedFilteredArrayElementStats.gc_depth,
            ],
            outdir = stats_out_dir + "/array_elements_for_transcriptome_discovery/",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_154_FinalizeAlignedArrayElementStats {
        input:
            files = [
                t_116_AlignedArrayElementStats.raw_stats,
                t_116_AlignedArrayElementStats.summary_stats,
                t_116_AlignedArrayElementStats.first_frag_qual,
                t_116_AlignedArrayElementStats.last_frag_qual,
                t_116_AlignedArrayElementStats.first_frag_gc_content,
                t_116_AlignedArrayElementStats.last_frag_gc_content,
                t_116_AlignedArrayElementStats.acgt_content_per_cycle,
                t_116_AlignedArrayElementStats.insert_size,
                t_116_AlignedArrayElementStats.read_length_dist,
                t_116_AlignedArrayElementStats.indel_distribution,
                t_116_AlignedArrayElementStats.indels_per_cycle,
                t_116_AlignedArrayElementStats.coverage_distribution,
                t_116_AlignedArrayElementStats.gc_depth,
            ],
            outdir = stats_out_dir + "/aligned_array_elements/",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_155_FinalizeBaselineArrayElementStats {
        input:
            files = [
                t_115_BaselineArrayElementStats.raw_stats,
                t_115_BaselineArrayElementStats.summary_stats,
                t_115_BaselineArrayElementStats.first_frag_qual,
                t_115_BaselineArrayElementStats.last_frag_qual,
                t_115_BaselineArrayElementStats.first_frag_gc_content,
                t_115_BaselineArrayElementStats.last_frag_gc_content,
                t_115_BaselineArrayElementStats.acgt_content_per_cycle,
                t_115_BaselineArrayElementStats.insert_size,
                t_115_BaselineArrayElementStats.read_length_dist,
                t_115_BaselineArrayElementStats.indel_distribution,
                t_115_BaselineArrayElementStats.indels_per_cycle,
                t_115_BaselineArrayElementStats.coverage_distribution,
                t_115_BaselineArrayElementStats.gc_depth,
            ],
            outdir = stats_out_dir + "/baseline_array_elements/",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_156_FinalizeCCSLongbowStats {
        input:
            files = [
                t_120_CCS_longbow_stats.summary_stats,
                t_120_CCS_longbow_stats.array_length_counts_plot_png,
                t_120_CCS_longbow_stats.array_length_counts_plot_svg,
                t_120_CCS_longbow_stats.ligation_heatmap_nn_png,
                t_120_CCS_longbow_stats.ligation_heatmap_nn_svg,
                t_120_CCS_longbow_stats.ligation_heatmap_png,
                t_120_CCS_longbow_stats.ligation_heatmap_svg,
                t_120_CCS_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_120_CCS_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_120_CCS_longbow_stats.ligation_heatmap_reduced_png,
                t_120_CCS_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/CCS_Corrected/",
            keyfile = keyfile
    }
    call FF.FinalizeToDir as t_157_FinalizeReclaimableLongbowStats {
        input:
            files = [
                t_121_Reclaimable_longbow_stats.summary_stats,
                t_121_Reclaimable_longbow_stats.array_length_counts_plot_png,
                t_121_Reclaimable_longbow_stats.array_length_counts_plot_svg,
                t_121_Reclaimable_longbow_stats.ligation_heatmap_nn_png,
                t_121_Reclaimable_longbow_stats.ligation_heatmap_nn_svg,
                t_121_Reclaimable_longbow_stats.ligation_heatmap_png,
                t_121_Reclaimable_longbow_stats.ligation_heatmap_svg,
                t_121_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_121_Reclaimable_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_121_Reclaimable_longbow_stats.ligation_heatmap_reduced_png,
                t_121_Reclaimable_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/CCS_Reclaimable/",
            keyfile = keyfile
    }
    call FF.FinalizeToDir as t_158_FinalizeReclaimedLongbowStats {
        input:
            files = [
                t_122_Reclaimed_longbow_stats.summary_stats,
                t_122_Reclaimed_longbow_stats.array_length_counts_plot_png,
                t_122_Reclaimed_longbow_stats.array_length_counts_plot_svg,
                t_122_Reclaimed_longbow_stats.ligation_heatmap_nn_png,
                t_122_Reclaimed_longbow_stats.ligation_heatmap_nn_svg,
                t_122_Reclaimed_longbow_stats.ligation_heatmap_png,
                t_122_Reclaimed_longbow_stats.ligation_heatmap_svg,
                t_122_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_122_Reclaimed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_122_Reclaimed_longbow_stats.ligation_heatmap_reduced_png,
                t_122_Reclaimed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/CCS_Reclaimed/",
            keyfile = keyfile
    }
    call FF.FinalizeToDir as t_159_FinalizeOverallLongbowStats {
        input:
            files = [
                t_125_Overall_longbow_stats.summary_stats,
                t_125_Overall_longbow_stats.array_length_counts_plot_png,
                t_125_Overall_longbow_stats.array_length_counts_plot_svg,
                t_125_Overall_longbow_stats.ligation_heatmap_nn_png,
                t_125_Overall_longbow_stats.ligation_heatmap_nn_svg,
                t_125_Overall_longbow_stats.ligation_heatmap_png,
                t_125_Overall_longbow_stats.ligation_heatmap_svg,
                t_125_Overall_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_125_Overall_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_125_Overall_longbow_stats.ligation_heatmap_reduced_png,
                t_125_Overall_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/Overall/",
            keyfile = keyfile
    }
    call FF.FinalizeToDir as t_160_FinalizeAllPassedLongbowStats {
        input:
            files = [
                t_123_Passed_longbow_stats.summary_stats,
                t_123_Passed_longbow_stats.array_length_counts_plot_png,
                t_123_Passed_longbow_stats.array_length_counts_plot_svg,
                t_123_Passed_longbow_stats.ligation_heatmap_nn_png,
                t_123_Passed_longbow_stats.ligation_heatmap_nn_svg,
                t_123_Passed_longbow_stats.ligation_heatmap_png,
                t_123_Passed_longbow_stats.ligation_heatmap_svg,
                t_123_Passed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_123_Passed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_123_Passed_longbow_stats.ligation_heatmap_reduced_png,
                t_123_Passed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/All_Longbow_Passed/",
            keyfile = keyfile
    }
    call FF.FinalizeToDir as t_161_FinalizeAllPassedLongbowStats {
        input:
            files = [
                t_124_Failed_longbow_stats.summary_stats,
                t_124_Failed_longbow_stats.array_length_counts_plot_png,
                t_124_Failed_longbow_stats.array_length_counts_plot_svg,
                t_124_Failed_longbow_stats.ligation_heatmap_nn_png,
                t_124_Failed_longbow_stats.ligation_heatmap_nn_svg,
                t_124_Failed_longbow_stats.ligation_heatmap_png,
                t_124_Failed_longbow_stats.ligation_heatmap_svg,
                t_124_Failed_longbow_stats.ligation_heatmap_nn_reduced_png,
                t_124_Failed_longbow_stats.ligation_heatmap_nn_reduced_svg,
                t_124_Failed_longbow_stats.ligation_heatmap_reduced_png,
                t_124_Failed_longbow_stats.ligation_heatmap_reduced_svg,
            ],
            outdir = stats_out_dir + "/longbow_stats/All_Longbow_Failed/",
            keyfile = keyfile
    }

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good:
    call FF.WriteCompletionFile as t_162_WriteCompletionFile {
        input:
            outdir = base_out_dir + "/",
            keyfile = keyfile
    }
}
