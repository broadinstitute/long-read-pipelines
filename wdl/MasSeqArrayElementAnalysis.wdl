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

workflow MasSeqArrayElementAnalysis {

    meta {
        description : "This workflow is designed to process already segmented and extracted array elements from the MASSeq v2 protocol and analyze them in the same manner as the main workflow."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {

        File array_element_bam
        String sample_name

        File transcriptome_quant_fasta
        File transcriptome_quant_fasta_index
        File transcriptome_quant_fasta_dict
        File transcriptome_quant_annotation_gtf

        Boolean force_anndata_gencode_overwrite = false

        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/MasSeqArrayElementAnalysis"

        File ref_fasta =  "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
        File ref_fasta_index = "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai"
        File ref_fasta_dict = "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict"
        File genome_annotation_gtf = "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.primary_assembly.annotation.gtf"

        File transcriptome_ref_fasta =  "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.pc_transcripts.fa"
        File transcriptome_ref_fasta_index = "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.pc_transcripts.fa.fai"
        File transcriptome_ref_fasta_dict = "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.pc_transcripts.dict"

        File intervals_of_interest = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/gencode.37.TCR_intervals.tsv"
        String interval_overlap_name = "is_tcr_overlapping"
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

    call PB.PBIndex as t_02_IndexArrayElements {
        input:
            bam = array_element_bam
    }

    call PB.ShardLongReads as t_03_ShardArrayElements {
        input:
            unaligned_bam = array_element_bam,
            unaligned_pbi = t_02_IndexArrayElements.pbindex,
            prefix = sample_name + "_ArrayElements_shard",
            num_shards = 300,
            runtime_attr_override = fast_network_attrs_preemptible,
    }

    # To properly count our transcripts we must throw away the non-primary and unaligned reads:
    RuntimeAttr filterReadsAttrs = object {
        cpu_cores: 4,
        preemptible_tries: 0
    }

    scatter (sharded_array_elements in t_03_ShardArrayElements.unmapped_shards) {

        ###############################################################
        # Genome alignments:
        call AR.Minimap2 as t_04_AlignArrayElementsToGenome {
            input:
                reads      = [ sharded_array_elements ],
                ref_fasta  = ref_fasta,
                map_preset = "splice:hq"
        }

        # We need to restore the annotations we created with the 10x tool to the aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_05_RestoreAnnotationsToGenomeAlignedBam {
            input:
                annotated_bam_file = sharded_array_elements,
                aligned_bam_file = t_04_AlignArrayElementsToGenome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        ###############################################################
        # Transcriptome alignments:
        call AR.Minimap2 as t_06_AlignArrayElementsToTranscriptome {
            input:
                reads      = [ sharded_array_elements ],
                ref_fasta  = transcriptome_quant_fasta,
                map_preset = "asm20"
        }
        # We need to restore the annotations we created with the 10x tool to the aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_07_RestoreAnnotationsToTranscriptomeAlignedBam {
            input:
                annotated_bam_file = sharded_array_elements,
                aligned_bam_file = t_06_AlignArrayElementsToTranscriptome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        call Utils.FilterReadsBySamFlags as t_08_RemoveUnmappedAndNonPrimaryReads {
            input:
                bam = t_07_RestoreAnnotationsToTranscriptomeAlignedBam.output_bam,
                sam_flags = "2308",
                prefix = sample_name + "_ArrayElements_Annotated_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }
        # Filter reads with no UMI tag:
        call Utils.FilterReadsWithTagValues as t_09_FilterReadsWithNoUMI {
            input:
                bam = t_08_RemoveUnmappedAndNonPrimaryReads.output_bam,
                tag = "ZU",
                value_to_remove = ".",
                prefix = sample_name + "_ArrayElements_Annotated_Aligned_PrimaryOnly_WithUMIs",
                runtime_attr_override = filterReadsAttrs
        }
        # Copy the contig to a tag.
        # By this point in the pipeline, array elements are aligned to a transcriptome, so this tag will
        # actually indicate the transcript to which each array element aligns.
        call TENX.CopyContigNameToReadTag as t_10_CopyContigNameToReadTag {
            input:
                aligned_bam_file = t_09_FilterReadsWithNoUMI.output_bam,
                prefix = sample_name + "_ArrayElements_Annotated_Aligned_PrimaryOnly_WithUMIs"
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

    # Merge all CCS bams together for this Subread BAM:
    RuntimeAttr merge_extra_cpu_attrs = object {
        cpu_cores: 4
    }
    call Utils.MergeBams as t_11_MergeGenomeAlignedExtractedArrayElements { input: bams = t_05_RestoreAnnotationsToGenomeAlignedBam.output_bam, prefix = sample_name + "_array_elements_longbow_extracted_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }

    # Now we merge together our TX-ome aligned stuff:
    call Utils.MergeBams as t_12_MergeTranscriptomeAlignedExtractedArrayElements { input: bams = t_07_RestoreAnnotationsToTranscriptomeAlignedBam.output_bam, prefix = sample_name + "_array_elements_longbow_extracted_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_13_MergePrimaryTranscriptomeAlignedArrayElements { input: bams = t_10_CopyContigNameToReadTag.output_bam, prefix = sample_name + "_array_elements_longbow_extracted_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    ##########
    # Quantify Transcripts:
    ##########

    call UMI_TOOLS.Run_Group as t_14_UMIToolsGroup {
        input:
            aligned_transcriptome_reads = t_13_MergePrimaryTranscriptomeAlignedArrayElements.merged_bam,
            aligned_transcriptome_reads_index = t_13_MergePrimaryTranscriptomeAlignedArrayElements.merged_bai,
            do_per_cell = true,
            prefix = "~{sample_name}_umi_tools_group"
    }

    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_15_CreateCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_14_UMIToolsGroup.output_bam,
            prefix = "~{sample_name}_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromTsv as t_16_CreateCountMatrixAnndataFromTsv {
        input:
            count_matrix_tsv = t_15_CreateCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = transcriptome_quant_annotation_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            force_anndata_gencode_overwrite = force_anndata_gencode_overwrite,
            prefix = "~{sample_name}_gene_tx_expression_count_matrix"
    }

    ######################################################################
    #             _____ _             _ _
    #            |  ___(_)_ __   __ _| (_)_______
    #            | |_  | | '_ \ / _` | | |_  / _ \
    #            |  _| | | | | | (_| | | |/ /  __/
    #            |_|   |_|_| |_|\__,_|_|_/___\___|
    #
    ######################################################################

    String base_out_dir = outdir + "/" + sample_name + "/" + t_01_WdlExecutionStartTimestamp.timestamp_string

    String array_element_dir = base_out_dir + "/aligned_array_elements"
    String quant_dir = base_out_dir + "/quant"


    ##############################################################################################################
    # Finalize the final annotated, aligned array elements:
    call FF.FinalizeToDir as t_17_FinalizeArrayElements {
        input:
            files = [
                t_11_MergeGenomeAlignedExtractedArrayElements.merged_bam,
                t_11_MergeGenomeAlignedExtractedArrayElements.merged_bai,
                t_12_MergeTranscriptomeAlignedExtractedArrayElements.merged_bam,
                t_12_MergeTranscriptomeAlignedExtractedArrayElements.merged_bai,
                t_13_MergePrimaryTranscriptomeAlignedArrayElements.merged_bam,
                t_13_MergePrimaryTranscriptomeAlignedArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = t_16_CreateCountMatrixAnndataFromTsv.transcript_gene_count_anndata_h5ad,
    }

    call FF.FinalizeToDir as t_18_FinalizeQuantInfo {
        input:
            files = [
                t_15_CreateCountMatrixFromAnnotatedBam.count_matrix,
                t_16_CreateCountMatrixAnndataFromTsv.transcript_gene_count_anndata_h5ad,
            ],
            outdir = quant_dir,
    }

    call FF.FinalizeToDir as t_19_FinalizeAnnDataPickles {
        input:
            files = t_16_CreateCountMatrixAnndataFromTsv.pickles,
            outdir = quant_dir,
            keyfile = t_16_CreateCountMatrixAnndataFromTsv.transcript_gene_count_anndata_h5ad,
    }

    call FF.FinalizeToDir as t_20_FinalizeUmiTools {
        input:
            files = [
                    t_14_UMIToolsGroup.output_bam,
                    t_14_UMIToolsGroup.output_tsv,
                    t_14_UMIToolsGroup.memory_log,
            ],
            outdir = quant_dir,
            keyfile = t_16_CreateCountMatrixAnndataFromTsv.transcript_gene_count_anndata_h5ad,
    }
}
