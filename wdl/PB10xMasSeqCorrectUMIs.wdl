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

workflow PB10xMasSeqCorrectUMIs {
    input {
        File eq_class_bam
        File st_gtf
        File combined_tx_eq_class_assignments
        File combined_tx_eq_class_defs
        File combined_gene_eq_class_assignments
        File combined_gene_eq_class_defs

        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/PB10xMasSeqCorrectUMIs"

        File genome_annotation_gtf = "gs://broad-dsde-methods-long-reads/resources/gencode_v37/gencode.v37.primary_assembly.annotation.gtf"

        File intervals_of_interest = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/gencode.37.TCR_intervals.tsv"

        Float min_read_quality = 0.0
        String interval_overlap_name = "is_tcr_overlapping"

        String SM
        String ID

        # Add a suffix here for our out directory so we can label runs:
        String out_dir_suffix = ""
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
    String DIR = SM + "." + ID

    call LONGBOW.Correct_UMI as t_103_LongbowCorrectUmi {
        input:
            bam = eq_class_bam,
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
            tx_equivalence_class_assignments = combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = SM + "_ccs_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_110_CreateCCSCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_109_CreateCCSCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = combined_gene_eq_class_defs,
            prefix = SM + "_ccs_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: 64}
    }

    # Create CLR count matrix and anndata:
    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_111_CreateCLRCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_106_GetCcsReclaimedReadsWithCorrectedUmis.bam_out,
            tx_equivalence_class_assignments = combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = SM + "_clr_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_112_CreateCLRCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_111_CreateCLRCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = combined_gene_eq_class_defs,
            prefix = SM + "_clr_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: 64}
    }

    # Create overall count matrix and anndata:
    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_113_CreateOverallCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = t_103_LongbowCorrectUmi.umi_corrected_bam,
            tx_equivalence_class_assignments = combined_tx_eq_class_assignments,
            umi_tag = "BX",
            prefix = SM + "_overall_gene_tx_expression_count_matrix"
    }

    call TX_POST.CreateCountMatrixAnndataFromEquivalenceClasses as t_114_CreateOverallCountMatrixAnndataFromEqClasses {
        input:
            count_matrix_tsv = t_113_CreateOverallCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = st_gtf,
            gencode_reference_gtf_file = genome_annotation_gtf,
            overlap_intervals = intervals_of_interest,
            overlap_interval_label = interval_overlap_name,
            tx_equivalence_class_assignments = combined_tx_eq_class_assignments,
            tx_equivalence_class_definitions = combined_tx_eq_class_defs,
            gene_equivalence_class_assignments = combined_gene_eq_class_assignments,
            gene_equivalence_class_definitions = combined_gene_eq_class_defs,
            prefix = SM + "_overall_gene_tx_expression_count_matrix",

            runtime_attr_override = object {mem_gb: 64}
    }


    call AM.SamtoolsStats as t_118_AlignedAnnotatedArrayElementsForQuantStats {
        input:
            bam = t_103_LongbowCorrectUmi.umi_corrected_bam
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
                combined_gene_eq_class_defs,
                combined_gene_eq_class_assignments,
                combined_tx_eq_class_defs,
                combined_tx_eq_class_assignments,
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


    call FF.FinalizeToDir as t_137_FinalizeAnnotatedArrayElements {
        input:
            files = [
                t_103_LongbowCorrectUmi.umi_corrected_bam,
                t_103_LongbowCorrectUmi.umi_corrected_bam_index,

                t_104_GetCcsCorrectedReadsWithCorrectedUmis.bam_out,
                t_106_GetCcsReclaimedReadsWithCorrectedUmis.bam_out
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

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good:
    call FF.WriteCompletionFile as t_162_WriteCompletionFile {
        input:
            outdir = base_out_dir + "/",
            keyfile = keyfile
    }
}
