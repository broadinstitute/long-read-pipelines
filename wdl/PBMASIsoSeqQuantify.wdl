version 1.0

######################################################################################
## A workflow that performs processing of MAS-ISO-seq data on a single sample from
## one or more flow cells. The workflow merges multiple samples into a single BAM
## prior to processing.
######################################################################################

import "tasks/Utils.wdl" as Utils
import "tasks/PBUtils.wdl" as PB
import "tasks/AlignReads.wdl" as AR
import "tasks/StringTie2.wdl" as ST2
import "tasks/Ten_X_Tool.wdl" as TENX
import "tasks/MASSeq.wdl" as MAS
import "tasks/TranscriptAnalysis/UMI_Tools.wdl" as UMI_TOOLS
import "tasks/TranscriptAnalysis/Preprocessing_Tasks.wdl" as TX_PRE
import "tasks/TranscriptAnalysis/Postprocessing_Tasks.wdl" as TX_POST
import "tasks/Finalize.wdl" as FF

workflow PBMASIsoSeqQuantify {
    meta {
        description : "Quantifies RNA isoform expression from given extracted, aligned, UMI- and CBC-annotated MAS-ISO-seq reads."
    }
    input {
        Array[File] aligned_ccs_bams
        Array[File] aligned_ccs_bais

        String participant_name

        File ref_map_file
        File ref_gtf

        File? intervals_of_interest
        String? interval_overlap_name

        String gcs_out_root_dir

        Boolean DEBUG_MODE = false
    }

    parameter_meta {
        aligned_ccs_bams: "GCS path to aligned MAS-ISO-seq CCS BAM files"
        aligned_ccs_bais: "GCS path to aligned MAS-ISO-seq CCS BAM file indices"
        participant_name: "name of the participant from whom these samples were obtained"

        ref_map_file:     "table indicating reference sequence and auxillary file locations"
        ref_gtf:          "GTF file to use for quantification"

        intervals_of_interest : "[optional] An interval list file containing intervals to mark in the final anndata object as overlapping the transcripts."
        interval_overlap_name : "[optional] The name of the annotation to add to the final anndata object for the column containing the overlap flag for transcripts that overlap intervals in the given intervals_of_interest file."

        gcs_out_root_dir: "GCS bucket to store the corrected/uncorrected reads, variants, and metrics files"

        DEBUG_MODE:         "[default valued] enables debugging tasks / subworkflows (default: false)"
    }

    Float min_ccs_read_quality = 0.0

    Map[String, String] ref_map = read_map(ref_map_file)

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_01_WdlExecutionStartTimestamp { input: }

    String outdir = if DEBUG_MODE then sub(gcs_out_root_dir, "/$", "") + "/PBMASIsoSeqQuantify/~{participant_name}/" + t_01_WdlExecutionStartTimestamp.timestamp_string else sub(gcs_out_root_dir, "/$", "") + "/PBMASIsoSeqQuantify/~{participant_name}/"


    # Gather across (potential multiple) input BAMs
    call Utils.MergeBams as t_02_MergeAllReads { input: bams = aligned_ccs_bams, prefix = participant_name }

    File bam = t_02_MergeAllReads.merged_bam
    File bai = t_02_MergeAllReads.merged_bai

#    call Utils.SubsetBam as SubsetToChr21 {
#        input:
#            bam = t_02_MergeAllReads.merged_bam,
#            bai = t_02_MergeAllReads.merged_bai,
#            locus = "chr21"
#    }
#    File bam = SubsetToChr21.subset_bam
#    File bai = SubsetToChr21.subset_bai

    # Restore the tags we have renamed in PBFlowcell:
    call MAS.RestoreSingleCellBamTagsForMasIsoSeqV0 as t_03_RestoreSingleCellBamTagsForMasIsoSeqV0 {
        input:
            bam = bam,
            prefix = participant_name + ".tag_restored"
    }

    # Remove unmapped, secondary, supplementary, mq0, length > 15kb, end clips > 1kb
    call MAS.FilterMasSeqReads as t_03_AlignmentFilterArrayElements {
        input:
            input_bam = t_03_RestoreSingleCellBamTagsForMasIsoSeqV0.bam_out,
            input_bai = t_03_RestoreSingleCellBamTagsForMasIsoSeqV0.bai,
            prefix = participant_name + ".alignment_filtered",
    }

    call PB.PBIndex as t_04_IndexFilteredReads { input: bam = t_03_AlignmentFilterArrayElements.bam }
    File pbi = t_04_IndexFilteredReads.pbi

    # Create a novel transcriptome:
    call ST2.Quantify as t_05_Stringtie2_Quantify {
        input:
            aligned_bam = t_03_AlignmentFilterArrayElements.bam,
            aligned_bai = t_03_AlignmentFilterArrayElements.bai,
            gtf = ref_gtf,
            keep_retained_introns = false,
            prefix = participant_name + "_StringTie2_Quantify",
    }

    call ST2.ExtractTranscriptSequences as t_06_Stringtie2_ExtractTranscriptSequences {
        input:
            ref_fasta = ref_map['fasta'],
            ref_fasta_fai = ref_map['fai'],
            gtf = t_05_Stringtie2_Quantify.st_gtf,
            prefix = participant_name + "_StringTie2_ExtractTranscriptSequences",
    }

    call ST2.CompareTranscriptomes as t_07_Stringtie2_CompareTranscriptomes {
        input:
            guide_gtf = ref_gtf,
            new_gtf = t_05_Stringtie2_Quantify.st_gtf,
            prefix = participant_name + "_StringTie2_CompareTranscriptome",
    }

    # Now create gff comparisons for the MEP method:
    call TX_PRE.GffCompare as t_08_GffCompareStringtie2toGencode {
        input:
            gff_ref = t_05_Stringtie2_Quantify.st_gtf,
            gff_query = ref_gtf,
            ref_fasta = ref_map['fasta'],
            ref_fasta_index = ref_map['fai'],
            prefix = participant_name + "_st2_vs_gencode"
    }
    call TX_PRE.GffCompare as t_09_GffCompareGencodetoStringtie2 {
        input:
            gff_ref = ref_gtf,
            gff_query = t_05_Stringtie2_Quantify.st_gtf,
            ref_fasta = ref_map['fasta'],
            ref_fasta_index = ref_map['fai'],
            prefix = participant_name + "_gencode_vs_st2"
    }

    # Restore original read names to reads to ease downstream processing:
    call TX_PRE.RestoreOriginalReadNames as t_10_RestoreOriginalReadNames {
        input:
            bam = t_03_AlignmentFilterArrayElements.bam,
            prefix =  participant_name + ".array_elements.alignment_filtered.original_read_names",
    }

    # Break Our BAM into fixed number of shards
    # TODO: RE-enable filtering for v2 of the pipeline:
#    Array[String] default_filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']
    call Utils.MakeChrIntervalList as t_11_MakeChrIntervalList { input: ref_dict = ref_map['dict'], filter = [] }

    # Now we have to align the array elements to the new transcriptome.
    scatter (c in t_11_MakeChrIntervalList.chrs) {
        String contig = c[0]

        call Utils.SubsetBam as t_12_SubsetBam { input: bam = t_10_RestoreOriginalReadNames.bam_out, bai = t_10_RestoreOriginalReadNames.bai_out, locus = contig }

        # Create a GFF file:
        call TX_PRE.ConvertSplicedBamToGff as t_13_ConvertSplicedBamToGff {
            input:
                bam = t_12_SubsetBam.subset_bam
        }

        # Compare GFF files:
        call TX_PRE.GffCompare as t_14_GffCompareStringtie2toMasSeqReads {
            input:
                gff_ref = t_05_Stringtie2_Quantify.st_gtf,
                gff_query = t_13_ConvertSplicedBamToGff.gff,
                ref_fasta = ref_map['fasta'],
                ref_fasta_index = ref_map['fai'],
                prefix = participant_name + "_" + contig + "_st2_vs_reads"
        }

        call TX_PRE.GffCompare as t_15_GffCompareGencodetoMasSeqReads {
            input:
                gff_ref = ref_gtf,
                gff_query = t_13_ConvertSplicedBamToGff.gff,
                ref_fasta = ref_map['fasta'],
                ref_fasta_index = ref_map['fai'],
                prefix = participant_name + "_" + contig + "_gencode_vs_reads"
        }

        # Create the comparison graph and tsv files:
        call TX_POST.QuantifyGffComparison as t_16_QuantifyGffComparison {
            input:
                genome_gtf = ref_gtf,
                st2_gencode_refmap = t_08_GffCompareStringtie2toGencode.refmap,
                st2_gencode_tmap = t_08_GffCompareStringtie2toGencode.tmap,
                st2_read_refmap = t_14_GffCompareStringtie2toMasSeqReads.refmap,
                st2_read_tmap = t_14_GffCompareStringtie2toMasSeqReads.tmap,
                gencode_st2_refmap = t_09_GffCompareGencodetoStringtie2.refmap,
                gencode_st2_tmap = t_09_GffCompareGencodetoStringtie2.tmap,
                gencode_read_refmap = t_15_GffCompareGencodetoMasSeqReads.refmap,
                gencode_read_tmap = t_15_GffCompareGencodetoMasSeqReads.tmap,
                prefix = participant_name + ".array_elements_" + contig
        }
    }

    # Merge our tx equivalance classes assignments and eq classes:
    call TX_POST.CombineEqClassFiles as t_17_CombineEqClassFiles {
        input:
            gene_eq_class_definitions = t_16_QuantifyGffComparison.gene_eq_class_labels_file,
            gene_assignment_files = t_16_QuantifyGffComparison.gene_assignments_file,
            equivalence_class_definitions = t_16_QuantifyGffComparison.tx_equivalence_class_labels_file,
            equivalence_classes = t_16_QuantifyGffComparison.tx_equivalence_class_file,
            prefix = participant_name + ".array_elements"
    }

    # Copy the eq class info to our annotated files so we can run UMI Cover:
    call TX_POST.CopyEqClassInfoToTag as t_18_CopyEqClassInfoToTag {
        input:
            bam = t_10_RestoreOriginalReadNames.bam_out,
            eq_class_file = t_17_CombineEqClassFiles.combined_tx_eq_class_assignments,
            prefix = participant_name + ".array_elements.alignment_filtered.eq_class_tags"
    }

    # Now we go in and fix our UMIs with set cover:
    call TX_PRE.CorrectUmisWithSetCover as t_19_CorrectUmisWithSetCover {
        input:
            bam = t_18_CopyEqClassInfoToTag.bam_out,
            prefix = participant_name + ".array_elements.alignment_filtered.for_quant"
    }

    # TODO: Debug memory requirement:
    Int anndata_mem_gb = 64

    # NOTE: We currently only allow CCS reads, so the task name here is correct.
    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_20_CreateCCSCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_19_CorrectUmisWithSetCover.corrected_umi_reads,
            tx_equivalence_class_assignments = t_17_CombineEqClassFiles.combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = participant_name + ".ccs_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_21_CreateCCSCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_20_CreateCCSCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = t_05_Stringtie2_Quantify.st_gtf,
            gencode_reference_gtf_file = ref_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = t_17_CombineEqClassFiles.combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = t_17_CombineEqClassFiles.combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = t_17_CombineEqClassFiles.combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = t_17_CombineEqClassFiles.combined_gene_eq_class_defs,
            prefix = participant_name + ".ccs_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: anndata_mem_gb}
    }

    ######################################################################
    #             _____ _             _ _
    #            |  ___(_)_ __   __ _| (_)_______
    #            | |_  | | '_ \ / _` | | |_  / _ \
    #            |  _| | | | | | (_| | | |/ /  __/
    #            |_|   |_|_| |_|\__,_|_|_/___\___|
    #
    ######################################################################

    # Key our finalizing off the file that will be created last so that we don't end up with partial runs:
    File keyfile = t_21_CreateCCSCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad

    String input_data_dir = outdir + "/combined_input_data"
    String adir = outdir + "/alignments"

    String quant_dir = outdir + "/quant"
    String meta_files_dir = outdir + "/meta_files"
    String intermediate_reads_dir = outdir + "/intermediate_reads/"

    call FF.FinalizeToFile as t_22_FinalizeBam { input: outdir = input_data_dir, file = bam, name = "~{participant_name}.bam", keyfile = keyfile }
    call FF.FinalizeToFile as t_23_FinalizeBai { input: outdir = input_data_dir, file = bai, name = "~{participant_name}.bam.bai", keyfile = keyfile }

    # TODO: Fix outputs / finalization!

    ##############################################################################################################
    # Finalize quantified expression:
    String eq_class_dir = quant_dir + "/eqivalence_classes"
    call FF.FinalizeToFile as t_24_FinalizeEqClassesGeneDefs { input: outdir = eq_class_dir, file = t_17_CombineEqClassFiles.combined_gene_eq_class_defs, keyfile = keyfile}
    call FF.FinalizeToFile as t_25_FinalizeEqClassesGeneAssignments { input: outdir = eq_class_dir, file = t_17_CombineEqClassFiles.combined_gene_eq_class_assignments, keyfile = keyfile}
    call FF.FinalizeToFile as t_26_FinalizeEqClassesTranscriptDefs { input: outdir = eq_class_dir, file = t_17_CombineEqClassFiles.combined_tx_eq_class_defs, keyfile = keyfile}
    call FF.FinalizeToFile as t_27_FinalizeEqClassesTranscriptAssignments { input: outdir = eq_class_dir, file = t_17_CombineEqClassFiles.combined_tx_eq_class_assignments, keyfile = keyfile}

    # CCS:
    call FF.FinalizeToFile as t_28_FinalizeCCSCountMatrix { input: outdir = quant_dir + "/CCS", file = t_20_CreateCCSCountMatrixFromAnnotatedBam.count_matrix, keyfile = keyfile}
    call FF.FinalizeToFile as t_29_FinalizeCCSCountMatrixAnnData { input: outdir = quant_dir + "/CCS", file = t_21_CreateCCSCountMatrixAnndataFromEqClasses.transcript_gene_count_anndata_h5ad, keyfile = keyfile}

    call FF.FinalizeToDir as t_30_FinalizeCCSRawQuantPickles {
        input:
            files = t_21_CreateCCSCountMatrixAnndataFromEqClasses.pickles,
            outdir = quant_dir + "/CCS",
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_31_FinalizeRefAndSt2Comparisons {
        input:
            files = [
                t_08_GffCompareStringtie2toGencode.refmap,
                t_08_GffCompareStringtie2toGencode.tmap,
                t_08_GffCompareStringtie2toGencode.tracking,
                t_08_GffCompareStringtie2toGencode.loci,
                t_08_GffCompareStringtie2toGencode.annotated_gtf,
                t_08_GffCompareStringtie2toGencode.stats,
                t_08_GffCompareStringtie2toGencode.log,

                t_09_GffCompareGencodetoStringtie2.refmap,
                t_09_GffCompareGencodetoStringtie2.tmap,
                t_09_GffCompareGencodetoStringtie2.tracking,
                t_09_GffCompareGencodetoStringtie2.loci,
                t_09_GffCompareGencodetoStringtie2.annotated_gtf,
                t_09_GffCompareGencodetoStringtie2.stats,
                t_09_GffCompareGencodetoStringtie2.log,
            ],
            outdir = quant_dir + "/gencode_and_stringtie2",
            keyfile = keyfile
    }

    # Finalize gene / tx assignment by contig:
    # NOTE: According to the scatter/gather documentation in the WDL spec, this will work correctly
    #       (https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#scatter--gather)
    scatter (i in range(length(t_11_MakeChrIntervalList.chrs))) {

        # NOTE: This must be `the_contg` because all variable names in WDL must be unique!
        String the_contig = t_11_MakeChrIntervalList.chrs[i][0]

        call FF.FinalizeToDir as t_32_FinalizeTxAndGeneAssignmentsByContig {
            input:
                files = [
                    t_14_GffCompareStringtie2toMasSeqReads.refmap[i],
                    t_14_GffCompareStringtie2toMasSeqReads.tmap[i],
                    t_14_GffCompareStringtie2toMasSeqReads.tracking[i],
                    t_14_GffCompareStringtie2toMasSeqReads.loci[i],
                    t_14_GffCompareStringtie2toMasSeqReads.annotated_gtf[i],
                    t_14_GffCompareStringtie2toMasSeqReads.stats[i],
                    t_14_GffCompareStringtie2toMasSeqReads.log[i],

                    t_15_GffCompareGencodetoMasSeqReads.refmap[i],
                    t_15_GffCompareGencodetoMasSeqReads.tmap[i],
                    t_15_GffCompareGencodetoMasSeqReads.tracking[i],
                    t_15_GffCompareGencodetoMasSeqReads.loci[i],
                    t_15_GffCompareGencodetoMasSeqReads.annotated_gtf[i],
                    t_15_GffCompareGencodetoMasSeqReads.stats[i],
                    t_15_GffCompareGencodetoMasSeqReads.log[i],

                    t_16_QuantifyGffComparison.gene_assignments_file[i],
                    t_16_QuantifyGffComparison.gene_eq_class_labels_file[i],
                    t_16_QuantifyGffComparison.tx_equivalence_class_labels_file[i],
                    t_16_QuantifyGffComparison.tx_equivalence_class_file[i],
                    t_16_QuantifyGffComparison.graph_gpickle[i],
                ],
                outdir = quant_dir + "/by_contig/" + the_contig,
                keyfile = keyfile
        }
    }

    ##############################################################################################################
    # Finalize annotated, aligned array elements:

    call FF.FinalizeToDir as t_33_FinalizeIntermediateAnnotatedArrayElements {
        input:
            files = [
                t_19_CorrectUmisWithSetCover.uncorrected_umi_reads
            ],
            outdir = intermediate_reads_dir,
            keyfile = keyfile
    }

    call FF.FinalizeToDir as t_34_FinalizeAnnotatedArrayElements {
        input:
            files = [
                t_19_CorrectUmisWithSetCover.corrected_umi_reads,
                t_19_CorrectUmisWithSetCover.corrected_umi_reads_index,
            ],
            outdir = adir,
            keyfile = keyfile
    }

    ##############################################################################################################
    # Finalize the discovered transcriptome:
    call FF.FinalizeToDir as t_35_FinalizeDiscoveredTranscriptome {
        input:
            files = [
                t_05_Stringtie2_Quantify.st_gtf,
                t_06_Stringtie2_ExtractTranscriptSequences.transcripts_fa,
                t_06_Stringtie2_ExtractTranscriptSequences.transcripts_fai,
                t_06_Stringtie2_ExtractTranscriptSequences.transcripts_dict,
                t_07_Stringtie2_CompareTranscriptomes.annotated_gtf,
                t_07_Stringtie2_CompareTranscriptomes.loci,
                t_07_Stringtie2_CompareTranscriptomes.stats,
                t_07_Stringtie2_CompareTranscriptomes.tracking,
                t_07_Stringtie2_CompareTranscriptomes.refmap,
                t_07_Stringtie2_CompareTranscriptomes.tmap,
            ],
            outdir = outdir + "/discovered_transcriptome",
            keyfile = keyfile
    }

    output {
        File merged_input_bam = t_22_FinalizeBam.gcs_path
        File merged_input_bai = t_23_FinalizeBai.gcs_path

        File gene_eq_class_defs = t_24_FinalizeEqClassesGeneDefs.gcs_path
        File gene_eq_class_assignments = t_25_FinalizeEqClassesGeneAssignments.gcs_path
        File transcript_eq_class_defs = t_26_FinalizeEqClassesTranscriptDefs.gcs_path
        File transcript_eq_class_assignments = t_27_FinalizeEqClassesTranscriptAssignments.gcs_path
        File raw_count_matrix = t_28_FinalizeCCSCountMatrix.gcs_path
        File count_matrix_anndata = t_29_FinalizeCCSCountMatrixAnnData.gcs_path

        String raw_quant_pickle_gcs_dir = t_30_FinalizeCCSRawQuantPickles.gcs_dir
        String ref_stringtie2_comparison_gcs_dir = t_31_FinalizeRefAndSt2Comparisons.gcs_dir
        Array[String] tx_gene_assignments_gcs_dir = t_32_FinalizeTxAndGeneAssignmentsByContig.gcs_dir
        String intermediate_read_gcs_dir = t_33_FinalizeIntermediateAnnotatedArrayElements.gcs_dir
        String annotated_array_element_gcs_dir = t_34_FinalizeAnnotatedArrayElements.gcs_dir
        String discovered_transcriptome_gcs_dir = t_35_FinalizeDiscoveredTranscriptome.gcs_dir
    }
}
