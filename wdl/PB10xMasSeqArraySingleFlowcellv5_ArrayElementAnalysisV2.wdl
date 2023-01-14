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

workflow PB10xMasSeqArraySingleFlowcellv5_ArrayElementAnalysisV2 {

    meta {
        description : "This workflow is designed to process data from the MASSeq v2 protocol and produce aligned reads that are ready for downstream analysis (e.g. transcript isoform identification).  It takes in a raw PacBio run folder location on GCS and produces a folder containing the aligned reads and other processed data."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String input_array_element_bam
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/PB10xMasSeqArraySingleFlowcellv5_ArrayElementAnalysisV2"

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

        String sample_name
        String sample_id
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

    ## No more preemption on this sharding - takes too long otherwise.
    RuntimeAttr disable_preemption_runtime_attrs = object {
        preemptible_tries: 0
    }

        # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_001_WdlExecutionStartTimestamp { input: }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    String SM  = sample_name
    String ID  = sample_id
    String DIR = SM + "." + ID

    String fbmrq_prefix = basename(input_array_element_bam, ".bam")

    Array[String] tags_to_preserve =  [ "CB", "JB", "JC", "JD", "JF", "JX", "RC", "RG", "SG", "XA", "XB", "XC", "XF", "XM", "XN", "XQ", "XU", "YC", "YG", "YK", "YN", "YP", "YQ", "YS", "YV", "ZS", "ZU", "ec", "fn", "ic", "im", "is", "it", "np", "pz", "rn", "rq", "sn", "we", "ws", "zm" ]

    call Utils.Bamtools as t_002_GetCcsReads {
        input:
            bamfile = input_array_element_bam,
            prefix = fbmrq_prefix + "_ccs_reads",
            cmd = "filter",
            args = '-tag "rq":">=' + min_read_quality + '"',
            runtime_attr_override = disable_preemption_runtime_attrs
    }

    call Utils.Bamtools as t_003_GetClrReads {
        input:
            bamfile = input_array_element_bam,
            prefix = fbmrq_prefix + "_clr_reads",
            cmd = "filter",
            args = '-tag "rq":"<' + min_read_quality + '"',
            runtime_attr_override = disable_preemption_runtime_attrs
    }

    ################################################################################################
    #    ____ ____ ____
    #   / ___/ ___/ ___|
    #  | |  | |   \___ \
    #  | |__| |___ ___) |
    #   \____\____|____/
    #
    ################################################################################################

    call PB.PBIndex as t_004_PbIndexCCSReads { input: bam = t_002_GetCcsReads.bam_out }
    call PB.ShardLongReads as t_005_ShardCcsReads {
        input:
            unaligned_bam = t_002_GetCcsReads.bam_out,
            unaligned_pbi = t_004_PbIndexCCSReads.pbindex,
            prefix = SM + "_ccs_shard",
            num_shards = 100,
    }
    scatter (ccs_shard_index in range(length(t_005_ShardCcsReads.unmapped_shards))) {
        File sharded_ccs_reads = t_005_ShardCcsReads.unmapped_shards[ccs_shard_index]

        # Now call the new Longbow Sift to remove individual reads that are incomplete / truncated / malformed:
        call LONGBOW.Sift as t_006_LongbowSiftCCSArrayElements {
            input:
                segmented_input_reads = sharded_ccs_reads,
                model = mas_seq_model,
                prefix = SM + "_ccs_array_elements_annotated_sifted_shard_" + ccs_shard_index,
                runtime_attr_override = object { cpu_cores: 4, mem_gb: 32, preemptible_tries: 0}
        }

        # Now that we've annotated the reads, we can pad the UMIs by a couple of bases to aid in the deduping:
        call LONGBOW.Pad as t_007_LongbowPadCCSArrayElementUMIs {
            input:
                reads = t_006_LongbowSiftCCSArrayElements.sifted_bam,
                model = mas_seq_model,
                tag_to_expand = "ZU",
                padding = ccs_umi_padding,
                prefix = SM + "_ccs_array_elements_annotated_umi_padded_shard_" + ccs_shard_index,
        }

        call LONGBOW.Pad as t_008_LongbowPadCCSArrayElementCBCs {
            input:
                reads = t_007_LongbowPadCCSArrayElementUMIs.padded_tag_bam,
                model = mas_seq_model,
                tag_to_expand = "CR",
                new_tag_dest = expanded_cbc_tag,
                padding = ccs_cbc_padding,
                prefix = SM + "_ccs_array_elements_annotated_cbc_padded_shard_" + ccs_shard_index,
        }


        # Now we should correct our barcodes based on the whitelist:
        call LONGBOW.Correct as t_009_LongbowCorrectCCSCorrectedArrayElementCBCs {
            input:
                reads = t_008_LongbowPadCCSArrayElementCBCs.padded_tag_bam,
                barcode_allow_list = cell_barcode_whitelist,
                model = mas_seq_model,
                ccs_lev_dist_threshold = ccs_lev_dist,
                clr_lev_dist_threshold = clr_lev_dist,
                prefix = SM + "_ccs_array_elements_annotated_padded_cbc_corrected_shard_" + ccs_shard_index,
                raw_barcode_tag = expanded_cbc_tag,
                corrected_barcode_tag = "CB",
        }

        call TENX.AdjustUmiSequenceWithAdapterAlignment as t_010_AdjustCCSUMIs {
            input:
                bam = t_009_LongbowCorrectCCSCorrectedArrayElementCBCs.corrected_barcodes_bam,
                short_read_umis = short_read_umis_tsv,
                prefix = SM + "_ccs_array_elements_annotated_padded_cbc_corrected_UMI_adjusted_shard_" + ccs_shard_index
        }

        call LONGBOW.Extract as t_011_LongbowExtractCcsArrayElements {
            input:
                bam = t_010_AdjustCCSUMIs.output_bam,
                prefix = SM + "_ccs_array_elements_cbc_umi_padded_extracted_shard_" + ccs_shard_index
        }

        call AR.Minimap2 as t_012_AlignCCSArrayElementsToGenome {
            input:
                reads      = [ t_011_LongbowExtractCcsArrayElements.extracted_bam ],
                ref_fasta  = ref_fasta,
                tags_to_preserve = tags_to_preserve,
                map_preset = "splice:hq",
                prefix = SM + "_ccs_array_elements_extracted_aligned_shard_" + ccs_shard_index,
                runtime_attr_override = object { mem_gb: 32 }
        }

        call LONGBOW.TagFix as t_013_LongbowTagfixAlignedCcsArrayElements {
            input:
                bam = t_012_AlignCCSArrayElementsToGenome.aligned_bam,
                prefix = SM + "_ccs_array_elements_extracted_aligned_tagfixed_shard_" + ccs_shard_index,
        }
    }

    call Utils.MergeBams as t_014_MergeCCSPLongbowTagfixedReads {
        input:
            bams = t_013_LongbowTagfixAlignedCcsArrayElements.tag_fixed_bam,
            prefix = SM + "_ccs_array_elements_extracted_aligned_tagfixed"
    }

    # Remove unmapped, secondary, supplementary, mq0, length > 15kb, end clips > 1kb
    call Utils.FilterMasSeqReadsWithGatk as t_015_AlignmentFilterForCcsArrayElements {
        input:
            bam_file = t_014_MergeCCSPLongbowTagfixedReads.merged_bam,
            bam_index = t_014_MergeCCSPLongbowTagfixedReads.merged_bai,
            prefix = SM + "_CCS_ArrayElements_Annotated_Aligned_PrimaryOnly",
            runtime_attr_override = disable_preemption_runtime_attrs
    }


    ################################################################################################
    #    ____ _     ____
    #   / ___| |   |  _ \
    #  | |   | |   | |_) |
    #  | |___| |___|  _ <
    #   \____|_____|_| \_\
    #
    ################################################################################################

    call PB.PBIndex as t_016_PbIndexCLRReads { input: bam = t_003_GetClrReads.bam_out }
    call PB.ShardLongReads as t_017_ShardClrReads {
        input:
            unaligned_bam = t_003_GetClrReads.bam_out,
            unaligned_pbi = t_016_PbIndexCLRReads.pbindex,
            prefix = SM + "_clr_shard",
            num_shards = 100,
    }
    scatter (clr_shard_index in range(length(t_017_ShardClrReads.unmapped_shards))) {
        File sharded_clr_reads = t_017_ShardClrReads.unmapped_shards[clr_shard_index]

        # Now call the new Longbow Sift to remove individual reads that are incomplete / truncated / malformed:
        call LONGBOW.Sift as t_018_LongbowSiftClrArrayElements {
            input:
                segmented_input_reads = sharded_clr_reads,
                model = mas_seq_model,
                prefix = SM + "_clr_array_elements_annotated_sifted_shard_" + clr_shard_index,
                runtime_attr_override = object { cpu_cores: 4, mem_gb: 32, preemptible_tries: 0}
        }

        # Now that we've annotated the reads, we can pad the UMIs by a couple of bases to aid in the deduping:
        call LONGBOW.Pad as t_019_LongbowPadClrArrayElementUMIs {
            input:
                reads = t_018_LongbowSiftClrArrayElements.sifted_bam,
                model = mas_seq_model,
                tag_to_expand = "ZU",
                padding = clr_umi_padding,
                prefix = SM + "_clr_array_elements_annotated_umi_padded_shard_" + clr_shard_index,
        }

        call LONGBOW.Pad as t_020_LongbowPadClrArrayElementCBCs {
            input:
                reads = t_019_LongbowPadClrArrayElementUMIs.padded_tag_bam,
                model = mas_seq_model,
                tag_to_expand = "CR",
                new_tag_dest = expanded_cbc_tag,
                padding = ccs_cbc_padding,
                prefix = SM + "_clr_array_elements_annotated_cbc_padded_shard_" + clr_shard_index,
        }

        # Now we should correct our barcodes based on the whitelist:
        call LONGBOW.Correct as t_021_LongbowCorrectClrArrayElementCBCs {
            input:
                reads = t_020_LongbowPadClrArrayElementCBCs.padded_tag_bam,
                barcode_allow_list = cell_barcode_whitelist,
                model = mas_seq_model,
                ccs_lev_dist_threshold = ccs_lev_dist,
                clr_lev_dist_threshold = clr_lev_dist,
                prefix = SM + "_clr_array_elements_annotated_padded_cbc_corrected_shard_" + clr_shard_index,
                raw_barcode_tag = expanded_cbc_tag,
                corrected_barcode_tag = "CB",
        }

        call TENX.AdjustUmiSequenceWithAdapterAlignment as t_022_AdjustClrUMIs {
            input:
                bam = t_021_LongbowCorrectClrArrayElementCBCs.corrected_barcodes_bam,
                short_read_umis = short_read_umis_tsv,
                prefix = SM + "_clr_array_elements_annotated_padded_cbc_corrected_UMI_adjusted_shard_" + clr_shard_index
        }

        call LONGBOW.Extract as t_023_LongbowExtractClrArrayElements {
            input:
                bam = t_022_AdjustClrUMIs.output_bam,
                prefix = SM + "_clr_array_elements_cbc_umi_padded_extracted_shard_" + clr_shard_index,
        }

        call AR.Minimap2 as t_024_AlignClrArrayElementsToGenome {
            input:
                reads      = [ t_023_LongbowExtractClrArrayElements.extracted_bam ],
                ref_fasta  = ref_fasta,
                tags_to_preserve = tags_to_preserve,
                map_preset = "splice",
                prefix = SM + "_clr_reclaimed_array_elements_extracted_aligned_shard_" + clr_shard_index,
                runtime_attr_override = object { mem_gb: 32 }
        }

        call LONGBOW.TagFix as t_025_LongbowTagfixAlignedClrArrayElements {
            input:
                bam = t_024_AlignClrArrayElementsToGenome.aligned_bam,
                prefix = SM + "_clr_reclaimed_array_elements_extracted_aligned_tagfixed_shard_" + clr_shard_index,
        }
    }

    # Merge the CLR Reads:
    call Utils.MergeBams as t_026_MergeCLRTagFixedReads {
        input:
            bams = t_025_LongbowTagfixAlignedClrArrayElements.tag_fixed_bam,
            prefix = SM + "_clr_array_elements_extracted_aligned_tagfixed"
    }

    call Utils.FilterMasSeqReadsWithGatk as t_027_AlignmentFilterForReclaimedArrayElements {
        input:
            bam_file = t_026_MergeCLRTagFixedReads.merged_bam,
            bam_index = t_026_MergeCLRTagFixedReads.merged_bai,
            prefix = SM + "_Reclaimed_ArrayElements_Annotated_Aligned_PrimaryOnly",
            runtime_attr_override = disable_preemption_runtime_attrs
    }

    ################################################################################################
    #   _______  __                            ____  _
    #  |_   _\ \/ /      ___  _ __ ___   ___  |  _ \(_)___  ___ _____   _____ _ __ _   _
    #    | |  \  /_____ / _ \| '_ ` _ \ / _ \ | | | | / __|/ __/ _ \ \ / / _ \ '__| | | |
    #    | |  /  \_____| (_) | | | | | |  __/ | |_| | \__ \ (_| (_) \ V /  __/ |  | |_| |
    #    |_| /_/\_\     \___/|_| |_| |_|\___| |____/|_|___/\___\___/ \_/ \___|_|   \__, |
    #                                                                              |___/
    ################################################################################################

    # Merge all alignments together:
    call Utils.MergeBams as t_028_MergeAllAlignedAndFilteredArrayElements {
        input:
            bams = [t_015_AlignmentFilterForCcsArrayElements.bam, t_027_AlignmentFilterForReclaimedArrayElements.bam],
            prefix = SM + "_all_array_elements_aligned_for_txome_discovery"
    }

    call StringTie2.Quantify as t_029_ST2_Quant {
        input:
            aligned_bam = t_028_MergeAllAlignedAndFilteredArrayElements.merged_bam,
            aligned_bai = t_028_MergeAllAlignedAndFilteredArrayElements.merged_bai,
            gtf = genome_annotation_gtf,
            keep_retained_introns = false,
            prefix = SM + "_StringTie2_Quantify",
    }

    call StringTie2.ExtractTranscriptSequences as t_030_ST2_ExtractTranscriptSequences  {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_index,
            gtf = t_029_ST2_Quant.st_gtf,
            prefix = SM + "_StringTie2_ExtractTranscriptSequences",
    }

    call StringTie2.CompareTranscriptomes as t_031_ST2_CompareTranscriptomes {
        input:
            guide_gtf = genome_annotation_gtf,
            new_gtf = t_029_ST2_Quant.st_gtf,
            prefix = SM + "_StringTie2_CompareTranscriptome",
    }

    ################################################################################################
    #   _____ ___         ____ _                 _____      _                  _   _
    #  | ____/ _ \       / ___| | __ _ ___ ___  | ____|_  _| |_ _ __ __ _  ___| |_(_) ___  _ __
    #  |  _|| | | |_____| |   | |/ _` / __/ __| |  _| \ \/ / __| '__/ _` |/ __| __| |/ _ \| '_ \
    #  | |__| |_| |_____| |___| | (_| \__ \__ \ | |___ >  <| |_| | | (_| | (__| |_| | (_) | | | |
    #  |_____\__\_\      \____|_|\__,_|___/___/ |_____/_/\_\\__|_|  \__,_|\___|\__|_|\___/|_| |_|
    #
    ################################################################################################

    #################
    # Here we restore the original read names to the bam because we're hashing them with Longbow.segment:

    # Restore original read names to CCS reads:
    call TX_PRE.RestoreOriginalReadNames as t_032_RestoreCcsOriginalReadNames {
        input:
            bam = t_015_AlignmentFilterForCcsArrayElements.bam,
            prefix =  SM + "_CCS_cbc_annotated_array_elements_padded_original_names"
    }

    # Restore original read names to CLR reads:
    call TX_PRE.RestoreOriginalReadNames as t_033_RestoreClrOriginalReadNames {
        input:
            bam = t_027_AlignmentFilterForReclaimedArrayElements.bam,
            prefix =  SM + "_CLR_cbc_annotated_array_elements_padded_original_names"
    }

    # Merge Aligned CCS and Reclaimed reads together:
    call Utils.MergeBams as t_034_MergeAllAnnotatedArrayElementsWithOriginalNames {
        input:
            bams = [t_032_RestoreCcsOriginalReadNames.bam_out, t_033_RestoreClrOriginalReadNames.bam_out],
            prefix = SM + "_all_cbc_annotated_array_elements_padded_original_names"
    }

    #################
    # Now we have to split the reads again, process them into gff files, run gffcompare and then aggregate the results in a graph

    # We can actually compare the references without needing to scatter:
    call TX_PRE.GffCompare as t_035_GffCompareStringtie2toGencode {
        input:
            gff_ref = t_029_ST2_Quant.st_gtf,
            gff_query = genome_annotation_gtf,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
    }
    call TX_PRE.GffCompare as t_036_GffCompareGencodetoStringtie2 {
        input:
            gff_ref = genome_annotation_gtf,
            gff_query = t_029_ST2_Quant.st_gtf,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
    }

    # Split by contig:
    call TX_PRE.SplitBamByContig as t_037_SplitArrayElementsByContig {
        input:
            bam = t_034_MergeAllAnnotatedArrayElementsWithOriginalNames.merged_bam,
            prefix = SM + "_all_cbc_annotated_array_elements_padded_original_names"
    }

    # For each contig:
    scatter (i in range(length(t_037_SplitArrayElementsByContig.contig_bams))) {

        File contig_bam = t_037_SplitArrayElementsByContig.contig_bams[i]
        String contig_name = t_037_SplitArrayElementsByContig.contig_names[i]

        # Create a GFF file:
        call TX_PRE.ConvertSplicedBamToGff as t_038_ConvertSplicedBamToGff {
            input:
                bam = contig_bam
        }

        # Compare GFF files:
        call TX_PRE.GffCompare as t_039_GffCompareStringtie2toMasSeqReads {
            input:
                gff_ref = t_029_ST2_Quant.st_gtf,
                gff_query = t_038_ConvertSplicedBamToGff.gff,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
        }

        call TX_PRE.GffCompare as t_040_GffCompareGencodetoMasSeqReads {
            input:
                gff_ref = genome_annotation_gtf,
                gff_query = t_038_ConvertSplicedBamToGff.gff,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
        }

        # Create the comparison graph and tsv files:
        call TX_POST.QuantifyGffComparison as t_041_QuantifyGffComparison {
            input:
                genome_gtf = genome_annotation_gtf,
                st2_gencode_refmap = t_035_GffCompareStringtie2toGencode.refmap,
                st2_gencode_tmap = t_035_GffCompareStringtie2toGencode.tmap,
                st2_read_refmap = t_039_GffCompareStringtie2toMasSeqReads.refmap,
                st2_read_tmap = t_039_GffCompareStringtie2toMasSeqReads.tmap,
                gencode_st2_refmap = t_036_GffCompareGencodetoStringtie2.refmap,
                gencode_st2_tmap = t_036_GffCompareGencodetoStringtie2.tmap,
                gencode_read_refmap = t_040_GffCompareGencodetoMasSeqReads.refmap,
                gencode_read_tmap = t_040_GffCompareGencodetoMasSeqReads.tmap,
                prefix = SM + "_all_cbc_annotated_array_elements_padded_" + contig_name
        }
    }

    # Merge our tx equivalance classes assignments and eq classes:
    call TX_POST.CombineEqClassFiles as t_042_CombineEqClassFiles {
        input:
            gene_eq_class_definitions = t_041_QuantifyGffComparison.gene_eq_class_labels_file,
            gene_assignment_files = t_041_QuantifyGffComparison.gene_assignments_file,
            equivalence_class_definitions = t_041_QuantifyGffComparison.tx_equivalence_class_labels_file,
            equivalence_classes = t_041_QuantifyGffComparison.tx_equivalence_class_file,
            prefix = SM + "_all_cbc_annotated_array_elements_padded"
    }

    ################################################################################################
    #   ___            __                         ___                    _
    #  |_ _|___  ___  / _| ___  _ __ _ __ ___    / _ \ _   _  __ _ _ __ | |_
    #   | |/ __|/ _ \| |_ / _ \| '__| '_ ` _ \  | | | | | | |/ _` | '_ \| __|
    #   | |\__ \ (_) |  _| (_) | |  | | | | | | | |_| | |_| | (_| | | | | |_
    #  |___|___/\___/|_|  \___/|_|  |_| |_| |_|  \__\_\\__,_|\__,_|_| |_|\__|
    #
    ################################################################################################

    # Use old quant method here as a baseline for comparison:
    call TX_POST.CopyEqClassInfoToTag as t_043_CopyEqClassInfoToTag {
        input:
            bam = t_034_MergeAllAnnotatedArrayElementsWithOriginalNames.merged_bam,
            eq_class_file = t_042_CombineEqClassFiles.combined_tx_eq_class_assignments,
            prefix = SM + "_annotated_array_elements_for_quant_with_gene_names"
    }

    call LONGBOW.Correct_UMI as t_044_LongbowCorrectUmi {
        input:
            bam = t_043_CopyEqClassInfoToTag.bam_out,
            prefix = SM + "_annotated_array_elements_for_quant_with_gene_names"
    }

    # Because of how we're doing things, we need to pull out the CCS and CCS Reclaimed reads from the output of the
    # set cover correction:
    call Utils.Bamtools as t_045_GetCcsCorrectedReadsWithCorrectedUmis {
        input:
            bamfile = t_044_LongbowCorrectUmi.umi_corrected_bam,
            prefix = SM + "_annotated_array_elements_for_quant_with_gene_names.corrected_umis.CCS",
            cmd = "filter",
            args = '-tag "rq":">=' + min_read_quality + '"',
            runtime_attr_override = disable_preemption_runtime_attrs
    }
    call Utils.IndexBam as t_046_IndexCcsReadsWithCorrectedUmis {input: bam = t_045_GetCcsCorrectedReadsWithCorrectedUmis.bam_out }

    call Utils.Bamtools as t_047_GetCcsReclaimedReadsWithCorrectedUmis {
        input:
            bamfile = t_044_LongbowCorrectUmi.umi_corrected_bam,
            prefix = SM + "_annotated_array_elements_for_quant_with_gene_names.corrected_umis.CCS_Reclaimed",
            cmd = "filter",
            args = '-tag "rq":"<' + min_read_quality + '"',
            runtime_attr_override = disable_preemption_runtime_attrs
    }
    call Utils.IndexBam as t_048_IndexCcsReclaimedReadsWithCorrectedUmis {input: bam = t_047_GetCcsReclaimedReadsWithCorrectedUmis.bam_out }

    call UMI_TOOLS.Run_Group as t_049_UMIToolsGroup {
        input:
            aligned_transcriptome_reads = t_044_LongbowCorrectUmi.umi_corrected_bam,
            aligned_transcriptome_reads_index = t_044_LongbowCorrectUmi.umi_corrected_bam_index,
            do_per_cell = true,
            prefix = SM + "_annotated_array_elements_with_gene_names_with_umi_tools_group_correction"
    }

    # Create CCS count matrix and anndata:
    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_050_CreateCCSCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_045_GetCcsCorrectedReadsWithCorrectedUmis.bam_out,
            tx_equivalence_class_assignments = t_042_CombineEqClassFiles.combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = SM + "_ccs_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_051_CreateCCSCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_050_CreateCCSCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = t_029_ST2_Quant.st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = t_042_CombineEqClassFiles.combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = t_042_CombineEqClassFiles.combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = t_042_CombineEqClassFiles.combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = t_042_CombineEqClassFiles.combined_gene_eq_class_defs,
            prefix = SM + "_ccs_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: 64}
    }

    # Create CLR count matrix and anndata:
    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_052_CreateCLRCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_047_GetCcsReclaimedReadsWithCorrectedUmis.bam_out,
            tx_equivalence_class_assignments = t_042_CombineEqClassFiles.combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = SM + "_clr_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_053_CreateCLRCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_052_CreateCLRCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = t_029_ST2_Quant.st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = t_042_CombineEqClassFiles.combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = t_042_CombineEqClassFiles.combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = t_042_CombineEqClassFiles.combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = t_042_CombineEqClassFiles.combined_gene_eq_class_defs,
            prefix = SM + "_clr_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: 64}
    }

    # Create overall count matrix and anndata:
    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_054_CreateOverallCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_044_LongbowCorrectUmi.umi_corrected_bam,
            tx_equivalence_class_assignments = t_042_CombineEqClassFiles.combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = SM + "_overall_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_055_CreateOverallCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_054_CreateOverallCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = t_029_ST2_Quant.st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = t_042_CombineEqClassFiles.combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = t_042_CombineEqClassFiles.combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = t_042_CombineEqClassFiles.combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = t_042_CombineEqClassFiles.combined_gene_eq_class_defs,
            prefix = SM + "_overall_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: 64}
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

    # This seems to take longer to get to:
    File keyfile = t_049_UMIToolsGroup.output_tsv

    String base_out_dir = outdir + "/" + DIR + "/" + out_dir_suffix + "/" + t_001_WdlExecutionStartTimestamp.timestamp_string
    String stats_out_dir = base_out_dir + "/stats"
    String array_element_dir = base_out_dir + "/annotated_array_elements"
    String intermediate_reads_dir = base_out_dir + "/intermediate_reads"

    String meta_files_dir = base_out_dir + "/meta_files"

    String intermediate_array_reads_dir = intermediate_reads_dir + "/array_reads"
    String intermediate_array_elements_dir = intermediate_reads_dir + "/array_elements"

    String quant_dir = base_out_dir + "/quant"

    ##############################################################################################################
    # Finalize gene / tx assignments:
    call FF.FinalizeToDir as t_056_FinalizeEqClasses {
        input:
            files = [
                t_042_CombineEqClassFiles.combined_gene_eq_class_defs,
                t_042_CombineEqClassFiles.combined_gene_eq_class_assignments,
                t_042_CombineEqClassFiles.combined_tx_eq_class_defs,
                t_042_CombineEqClassFiles.combined_tx_eq_class_assignments,
            ],
            outdir = quant_dir + "/eqivalence_classes",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_057_FinalizeUmiToolsOutputs {
        input:
            files = [
                t_049_UMIToolsGroup.output_bam,
                t_049_UMIToolsGroup.output_tsv,
            ],
            outdir = quant_dir + "/UMITools",
            keyfile = keyfile
    }

    # CCS:
    call FF.FinalizeToDir as t_058_FinalizeCCSTxAndGeneAssignments {
        input:
            files = [
                t_050_CreateCCSCountMatrixFromAnnotatedBam.count_matrix,
                t_051_CreateCCSCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad,
            ],
            outdir = quant_dir + "/CCS",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_059_FinalizeCCSRawQuantPickles {
        input:
            files = t_051_CreateCCSCountMatrixAnndataFromEqClasses.pickles,
            outdir = quant_dir + "/CCS",
            keyfile = keyfile
    }

    # CLR:
    call FF.FinalizeToDir as t_060_FinalizeCLRTxAndGeneAssignments {
        input:
            files = [
                t_052_CreateCLRCountMatrixFromAnnotatedBam.count_matrix,
                t_053_CreateCLRCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad,
            ],
            outdir = quant_dir + "/CLR",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_061_FinalizeCLRRawQuantPickles {
        input:
            files = t_053_CreateCLRCountMatrixAnndataFromEqClasses.pickles,
            outdir = quant_dir + "/CLR",
            keyfile = keyfile
    }

    # Overall:
    call FF.FinalizeToDir as t_062_FinalizeOverallTxAndGeneAssignments {
        input:
            files = [
                t_054_CreateOverallCountMatrixFromAnnotatedBam.count_matrix,
                t_055_CreateOverallCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad,
            ],
            outdir = quant_dir + "/Overall",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_063_FinalizeOverallRawQuantPickles {
        input:
            files = t_055_CreateOverallCountMatrixAnndataFromEqClasses.pickles,
            outdir = quant_dir + "/Overall",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_064_FinalizeRefAndSt2Comparisons {
        input:
            files = [
                t_035_GffCompareStringtie2toGencode.refmap,
                t_035_GffCompareStringtie2toGencode.tmap,
                t_035_GffCompareStringtie2toGencode.tracking,
                t_035_GffCompareStringtie2toGencode.loci,
                t_035_GffCompareStringtie2toGencode.annotated_gtf,
                t_035_GffCompareStringtie2toGencode.stats,
                t_035_GffCompareStringtie2toGencode.log,

                t_036_GffCompareGencodetoStringtie2.refmap,
                t_036_GffCompareGencodetoStringtie2.tmap,
                t_036_GffCompareGencodetoStringtie2.tracking,
                t_036_GffCompareGencodetoStringtie2.loci,
                t_036_GffCompareGencodetoStringtie2.annotated_gtf,
                t_036_GffCompareGencodetoStringtie2.stats,
                t_036_GffCompareGencodetoStringtie2.log,
            ],
            outdir = quant_dir + "/gencode_and_stringtie2",
            keyfile = keyfile
    }

    # Finalize gene / tx assignment by contig:
    # NOTE: According to the scatter/gather documentation in the WDL spec, this will work correctly
    #       (https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#scatter--gather)
    scatter (i in range(length(t_037_SplitArrayElementsByContig.contig_bams))) {
        String contig = t_037_SplitArrayElementsByContig.contig_names[i]

        call FF.FinalizeToDir as t_065_FinalizeTxAndGeneAssignmentsByContig {
            input:
                files = [
                    t_039_GffCompareStringtie2toMasSeqReads.refmap[i],
                    t_039_GffCompareStringtie2toMasSeqReads.tmap[i],
                    t_039_GffCompareStringtie2toMasSeqReads.tracking[i],
                    t_039_GffCompareStringtie2toMasSeqReads.loci[i],
                    t_039_GffCompareStringtie2toMasSeqReads.annotated_gtf[i],
                    t_039_GffCompareStringtie2toMasSeqReads.stats[i],
                    t_039_GffCompareStringtie2toMasSeqReads.log[i],

                    t_040_GffCompareGencodetoMasSeqReads.refmap[i],
                    t_040_GffCompareGencodetoMasSeqReads.tmap[i],
                    t_040_GffCompareGencodetoMasSeqReads.tracking[i],
                    t_040_GffCompareGencodetoMasSeqReads.loci[i],
                    t_040_GffCompareGencodetoMasSeqReads.annotated_gtf[i],
                    t_040_GffCompareGencodetoMasSeqReads.stats[i],
                    t_040_GffCompareGencodetoMasSeqReads.log[i],

                    t_041_QuantifyGffComparison.gene_assignments_file[i],
                    t_041_QuantifyGffComparison.gene_eq_class_labels_file[i],
                    t_041_QuantifyGffComparison.tx_equivalence_class_labels_file[i],
                    t_041_QuantifyGffComparison.tx_equivalence_class_file[i],
                    t_041_QuantifyGffComparison.graph_gpickle[i],
                ],
                outdir = quant_dir + "/by_contig/" + contig,
                keyfile = keyfile
        }
    }

    call FF.FinalizeToDir as t_066_FinalizeAnnotatedArrayElements {
        input:
            files = [
                t_044_LongbowCorrectUmi.umi_corrected_bam,
                t_044_LongbowCorrectUmi.umi_corrected_bam_index,
                t_044_LongbowCorrectUmi.failed_umi_correction_bam,
                t_044_LongbowCorrectUmi.cached_read_loci,

                t_045_GetCcsCorrectedReadsWithCorrectedUmis.bam_out,
                t_047_GetCcsReclaimedReadsWithCorrectedUmis.bam_out,

                t_028_MergeAllAlignedAndFilteredArrayElements.merged_bam,
                t_028_MergeAllAlignedAndFilteredArrayElements.merged_bai
            ],
            outdir = array_element_dir,
            keyfile = keyfile
    }

    call FF.FinalizeToFile as t_067_FinalizeCcsArrayElementCorrectedUmiIndex {
        input:
            file = t_046_IndexCcsReadsWithCorrectedUmis.bai,
            outfile = array_element_dir + "/" + SM + "_annotated_array_elements_for_quant_with_gene_names.corrected_umis.CCS.bam.bai",
            keyfile = keyfile
    }

    call FF.FinalizeToFile as t_068_FinalizeCcsReclaimedArrayElementCorrectedUmiIndex {
        input:
            file = t_048_IndexCcsReclaimedReadsWithCorrectedUmis.bai,
            outfile = array_element_dir + "/" + SM + "_annotated_array_elements_for_quant_with_gene_names.corrected_umis.CCS_Reclaimed",
            keyfile = keyfile
    }

    scatter (ccsi2 in range(length(t_005_ShardCcsReads.unmapped_shards))) {

        call FF.FinalizeToFile as t_069_FinalizeCCSCBCcorrectionLogsToMeta {
            input:
                file = t_009_LongbowCorrectCCSCorrectedArrayElementCBCs.log[ccsi2],
                outfile = meta_files_dir + "/" + "ccs_cbc_correction_logs/" + SM + "_ccs_array_elements_annotated_padded_cbc_corrected_shard_" + ccsi2 + ".log",
                keyfile = keyfile
        }

        call FF.FinalizeToFile as t_070_FinalizeCCSUmiAdjustmentLogs {
            input:
                file = t_010_AdjustCCSUMIs.log[ccsi2],
                outfile = meta_files_dir + "/" + "umi_adjustment_logs_ccs/" + SM + "_ccs_array_elements_cbc_umi_padded_extracted_shard_" + ccsi2 + ".log",
                keyfile = keyfile
        }
    }

    scatter (clri2 in range(length(t_017_ShardClrReads.unmapped_shards))) {

        call FF.FinalizeToFile as t_071_FinalizeCLRCBCcorrectionLogsToMeta {
            input:
                file = t_021_LongbowCorrectClrArrayElementCBCs.log[clri2],
                outfile = meta_files_dir + "/" + "clr_cbc_correction_logs/" + SM + "_clr_array_elements_annotated_padded_cbc_corrected_shard_" + clri2 + ".log",
                keyfile = keyfile
        }

        call FF.FinalizeToFile as t_072_FinalizeCLRUmiAdjustmentLogs {
            input:
                file = t_022_AdjustClrUMIs.log[clri2],
                outfile = meta_files_dir + "/" + "umi_adjustment_logs_clr/" + SM + "_clr_array_elements_cbc_umi_padded_extracted_shard_" + clri2 + ".log",
                keyfile = keyfile
        }
    }

    ##############################################################################################################
    # Finalize the discovered transcriptome:
    if ( !is_SIRV_data ) {
        call FF.FinalizeToDir as t_073_FinalizeDiscoveredTranscriptome {
            input:
                files = [
                    t_029_ST2_Quant.st_gtf,
                    t_030_ST2_ExtractTranscriptSequences.transcripts_fa,
                    t_030_ST2_ExtractTranscriptSequences.transcripts_fai,
                    t_030_ST2_ExtractTranscriptSequences.transcripts_dict,
                    t_031_ST2_CompareTranscriptomes.annotated_gtf,
                    t_031_ST2_CompareTranscriptomes.loci,
                    t_031_ST2_CompareTranscriptomes.stats,
                    t_031_ST2_CompareTranscriptomes.tracking,
                    t_031_ST2_CompareTranscriptomes.refmap,
                    t_031_ST2_CompareTranscriptomes.tmap,
                ],
                outdir = base_out_dir + "/discovered_transcriptome",
                keyfile = keyfile
        }
    }
    # Finalize Longbow Sift stats:
    scatter (i_f1 in range(length(t_005_ShardCcsReads.unmapped_shards))) {
        call FF.FinalizeToFile as t_074_FinalizeCcsLongbowSiftStats {
            input:
                file = t_006_LongbowSiftCCSArrayElements.stats_tsv[i_f1],
                outfile = stats_out_dir + "/longbow_stats/sift/ccs/" + SM + "_ccs_array_elements_annotated_sifted_" + i_f1 + ".stats.tsv",
                keyfile = keyfile
        }
        call FF.FinalizeToFile as t_075_FinalizeCcsLongbowSiftSummaryStats {
            input:
                file = t_006_LongbowSiftCCSArrayElements.summary_stats_tsv[i_f1],
                outfile = stats_out_dir + "/longbow_stats/sift/ccs/" + SM + "_ccs_array_elements_annotated_sifted_" + i_f1 + ".summary_stats.tsv",
                keyfile = keyfile
        }
    }
    scatter (i_f2 in range(length(t_017_ShardClrReads.unmapped_shards))) {
        call FF.FinalizeToFile as t_076_FinalizeCcsReclaimedLongbowSiftStats {
            input:
                file = t_018_LongbowSiftClrArrayElements.stats_tsv[i_f2],
                outfile = stats_out_dir + "/longbow_stats/sift/ccs_reclaimed/" + SM + "_ccs_reclaimed_array_elements_annotated_sifted_" + i_f2 + ".stats.tsv",
                keyfile = keyfile
        }
        call FF.FinalizeToFile as t_077_FinalizeCcsReclaimedLongbowSiftSummaryStats {
            input:
                file = t_018_LongbowSiftClrArrayElements.summary_stats_tsv[i_f2],
                outfile = stats_out_dir + "/longbow_stats/sift/ccs_reclaimed/" + SM + "_ccs_reclaimed_array_elements_annotated_sifted_" + i_f2 + ".summary_stats.tsv",
                keyfile = keyfile
        }
    }

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good:
    call FF.WriteCompletionFile as t_078_WriteCompletionFile {
        input:
            outdir = base_out_dir + "/",
            keyfile = keyfile
    }
}
