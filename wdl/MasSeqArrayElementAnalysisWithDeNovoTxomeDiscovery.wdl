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

workflow MasSeqArrayElementAnalysisWithDeNovoTxomeDiscovery {

    meta {
        description : "This workflow is designed to process already segmented and extracted array elements from the MASSeq v2 protocol and analyze them in the same manner as the main workflow."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {

        File array_element_bam
        String sample_name

        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/MasSeqArrayElementAnalysisWithDeNovoTxomeDiscovery"

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
    }

    RuntimeAttr merge_extra_cpu_attrs = object {
        cpu_cores: 4
    }
    call Utils.MergeBams as t_06_MergeGenomeAlignedExtractedArrayElements { input: bams = t_05_RestoreAnnotationsToGenomeAlignedBam.output_bam, prefix = sample_name + "_array_elements_longbow_extracted_genome_aligned", runtime_attr_override = merge_extra_cpu_attrs }

    RuntimeAttr st2_big_mem = object {
        cpu_cores: 4,
        mem_gb: 16
    }
    call StringTie2.Quantify as t_07_ST2_Quant {
        input:
            aligned_bam = t_06_MergeGenomeAlignedExtractedArrayElements.merged_bam,
            aligned_bai = t_06_MergeGenomeAlignedExtractedArrayElements.merged_bai,
            gtf = genome_annotation_gtf,
            keep_retained_introns = false,
            prefix = sample_name + "_StringTie2_Quantify",
            runtime_attr_override = st2_big_mem,
    }

    call StringTie2.ExtractTranscriptSequences as t_08_ST2_ExtractTranscriptSequences  {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_index,
            gtf = t_07_ST2_Quant.st_gtf,
            prefix = sample_name + "_StringTie2_ExtractTranscriptSequences",
    }

    call StringTie2.CompareTranscriptomes as t_09_ST2_CompareTranscriptomes {
        input:
            guide_gtf = genome_annotation_gtf,
            new_gtf = t_07_ST2_Quant.st_gtf,
            prefix = sample_name + "_StringTie2_CompareTranscriptome",
    }

    scatter (sharded_array_elements in t_03_ShardArrayElements.unmapped_shards) {

        ###############################################################
        # Transcriptome alignments:
        call AR.Minimap2 as t_10_AlignArrayElementsToTranscriptome {
            input:
                reads      = [ sharded_array_elements ],
                ref_fasta  = t_08_ST2_ExtractTranscriptSequences.transcripts_fa,
                map_preset = "asm20"
        }
        # We need to restore the annotations we created with the 10x tool to the aligned reads.
        call TENX.RestoreAnnotationstoAlignedBam as t_11_RestoreAnnotationsToTranscriptomeAlignedBam {
            input:
                annotated_bam_file = sharded_array_elements,
                aligned_bam_file = t_10_AlignArrayElementsToTranscriptome.aligned_bam,
                tags_to_ignore = [],
                mem_gb = 8,
        }

        call Utils.FilterReadsBySamFlags as t_12_RemoveUnmappedAndNonPrimaryReads {
            input:
                bam = t_11_RestoreAnnotationsToTranscriptomeAlignedBam.output_bam,
                sam_flags = "2308",
                prefix = sample_name + "_ArrayElements_Annotated_Aligned_PrimaryOnly",
                runtime_attr_override = filterReadsAttrs
        }
        # Filter reads with no UMI tag:
        call Utils.FilterReadsWithTagValues as t_13_FilterReadsWithNoUMI {
            input:
                bam = t_12_RemoveUnmappedAndNonPrimaryReads.output_bam,
                tag = "ZU",
                value_to_remove = ".",
                prefix = sample_name + "_ArrayElements_Annotated_Aligned_PrimaryOnly_WithUMIs",
                runtime_attr_override = filterReadsAttrs
        }
        # Copy the contig to a tag.
        # By this point in the pipeline, array elements are aligned to a transcriptome, so this tag will
        # actually indicate the transcript to which each array element aligns.
        call TENX.CopyContigNameToReadTag as t_14_CopyContigNameToReadTag {
            input:
                aligned_bam_file = t_13_FilterReadsWithNoUMI.output_bam,
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

    # Now we merge together our TX-ome aligned stuff:
    call Utils.MergeBams as t_15_MergeTranscriptomeAlignedExtractedArrayElements { input: bams = t_11_RestoreAnnotationsToTranscriptomeAlignedBam.output_bam, prefix = sample_name + "_array_elements_longbow_extracted_tx_aligned", runtime_attr_override = merge_extra_cpu_attrs }
    call Utils.MergeBams as t_16_MergePrimaryTranscriptomeAlignedArrayElements { input: bams = t_14_CopyContigNameToReadTag.output_bam, prefix = sample_name + "_array_elements_longbow_extracted_tx_aligned_primary_alignments", runtime_attr_override = merge_extra_cpu_attrs }

    ##########
    # Quantify Transcripts:
    ##########

    call UMI_TOOLS.Run_Group as t_17_UMIToolsGroup {
        input:
            aligned_transcriptome_reads = t_16_MergePrimaryTranscriptomeAlignedArrayElements.merged_bam,
            aligned_transcriptome_reads_index = t_16_MergePrimaryTranscriptomeAlignedArrayElements.merged_bai,
            do_per_cell = true,
            prefix = "~{sample_name}_umi_tools_group"
    }

    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_18_CreateCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_17_UMIToolsGroup.output_bam,
            prefix = "~{sample_name}_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromTsv as t_19_CreateCountMatrixAnndataFromTsv {
        input:
            count_matrix_tsv = t_18_CreateCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = t_07_ST2_Quant.st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
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
    String transcriptome_dir = base_out_dir + "/discovered_transcriptome"

    ##############################################################################################################
    # Finalize the final annotated, aligned array elements:
    call FF.FinalizeToDir as t_20_FinalizeArrayElements {
        input:
            files = [
                t_06_MergeGenomeAlignedExtractedArrayElements.merged_bam,
                t_06_MergeGenomeAlignedExtractedArrayElements.merged_bai,
                t_15_MergeTranscriptomeAlignedExtractedArrayElements.merged_bam,
                t_15_MergeTranscriptomeAlignedExtractedArrayElements.merged_bai,
                t_16_MergePrimaryTranscriptomeAlignedArrayElements.merged_bam,
                t_16_MergePrimaryTranscriptomeAlignedArrayElements.merged_bai,
            ],
            outdir = array_element_dir,
            keyfile = t_19_CreateCountMatrixAnndataFromTsv.transcript_gene_count_anndata_h5ad,
    }

    call FF.FinalizeToDir as t_21_FinalizeQuantInfo {
        input:
            files = [
                t_18_CreateCountMatrixFromAnnotatedBam.count_matrix,
                t_19_CreateCountMatrixAnndataFromTsv.transcript_gene_count_anndata_h5ad,
            ],
            outdir = quant_dir,
    }

    call FF.FinalizeToDir as t_22_FinalizeAnnDataPickles {
        input:
            files = t_19_CreateCountMatrixAnndataFromTsv.pickles,
            outdir = quant_dir,
    }

    call FF.FinalizeToDir as t_23_FinalizeUmiTools {
        input:
            files = [
                    t_17_UMIToolsGroup.output_bam,
                    t_17_UMIToolsGroup.output_tsv,
                    t_17_UMIToolsGroup.memory_log,
            ],
            outdir = quant_dir,
            keyfile = t_19_CreateCountMatrixAnndataFromTsv.transcript_gene_count_anndata_h5ad,
    }

    call FF.FinalizeToDir as t_24_FinalizeDiscoveredTranscriptome {
        input:
            files = [
                t_07_ST2_Quant.st_gtf,
                t_08_ST2_ExtractTranscriptSequences.transcripts_fa,
                t_08_ST2_ExtractTranscriptSequences.transcripts_fai,
                t_08_ST2_ExtractTranscriptSequences.transcripts_dict,
                t_09_ST2_CompareTranscriptomes.annotated_gtf,
                t_09_ST2_CompareTranscriptomes.loci,
                t_09_ST2_CompareTranscriptomes.stats,
                t_09_ST2_CompareTranscriptomes.tracking,
                t_09_ST2_CompareTranscriptomes.refmap,
                t_09_ST2_CompareTranscriptomes.tmap,
            ],
            outdir = transcriptome_dir,
            keyfile = t_19_CreateCountMatrixAnndataFromTsv.transcript_gene_count_anndata_h5ad,
    }
}
