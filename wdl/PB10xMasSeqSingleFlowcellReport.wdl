version 1.0

import "tasks/JupyterNotebooks.wdl" as JUPYTER

workflow PB10xMasSeqSingleFlowcellReport {

    meta {
        description : "Create a report for a given MASSeq run which summarizes the results  using a given Jupyter Notebook template."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File notebook_template              = "gs://broad-dsde-methods-long-reads/covid-19-aziz/MAS-seq_QC_report_template-interactive.ipynb"

        File subreads_stats
        File ccs_reads_stats
        File array_elements_stats

        File ccs_report_file
        File ccs_bam_file

        File array_element_bam_file

        File ebr_element_marker_alignments
        File ebr_initial_section_alignments
        File ebr_final_section_alignments
        File ebr_bounds_file                = "gs://broad-dsde-methods-long-reads/covid-19-aziz/bounds_file_for_extraction.txt"

        File ten_x_metrics_file
        File rna_seq_metrics_file
        File workflow_dot_file              = "gs://broad-dsde-methods-long-reads/covid-19-aziz/PB10xMasSeqArraySingleFlowcell.dot"
    }

    parameter_meta {
        notebook_template : "Jupyter notebook MASSeq template to run with the given data to produce a MASSeq report."

        subreads_stats : "Samtools stats file created from the raw subreads from the PacBio instrument."
        ccs_reads_stats : "Samtools raw stats file created from the aligned CCS corrected reads from the PacBio instrument."
        array_elements_stats : "Samtools raw stats file created from the aligned MASSeq array elements."

        ccs_report_file : "CCS report file from the CCS run for the data from the PacBio instrument."
        ccs_bam_file : "Unaligned reads file in BAM format from the CCS process (pre-array splitting)."

        array_element_bam_file : "Aligned reads file in BAM format containing aligned MASSeq array elements as individual reads."

        ebr_element_marker_alignments : "Raw marker alignments file from ExtractBoundedReads for the data from this MASSeq run."
        ebr_initial_section_alignments : "Initial section alignments file from ExtractBoundedReads for the data from this MASSeq run."
        ebr_final_section_alignments : "Final section alignments file from ExtractBoundedReads for the data from this MASSeq run."
        ebr_bounds_file : "Text file containing two comma-separated known segment names on each line.  These entries define delimited sections that were extracted from the reads and treated as individual array elements."

        ten_x_metrics_file : "Stats file from the 10x tool run for the data in this MASSeq run."
        rna_seq_metrics_file : "Picard CollectRnaSeqMetrics metrics file created from the aligned MASSeq array elements."
        workflow_dot_file : "DOT file containing the representation of this workflow used to create and analyze the data.  This is included in the QC reports (the DOT file can be generated with womtool)."
    }

    ## NOTE: This assumes ONE file for both the raw input and the 10x array element stats!
    ##       This should be fixed in version 2.
    call JUPYTER.PB10xMasSeqSingleFlowcellReport as GenerateReport {
        input:
            notebook_template              = notebook_template,

            subreads_stats                 = subreads_stats,
            ccs_reads_stats                = ccs_reads_stats,
            array_elements_stats           = array_elements_stats,
            ccs_report_file                = ccs_report_file,

            ccs_bam_file                   = ccs_bam_file,
            array_element_bam_file         = array_element_bam_file,

            ebr_element_marker_alignments  = ebr_element_marker_alignments,
            ebr_initial_section_alignments = ebr_initial_section_alignments,
            ebr_final_section_alignments   = ebr_final_section_alignments,
            ebr_bounds_file                = ebr_bounds_file,

            ten_x_metrics_file             = ten_x_metrics_file,
            rna_seq_metrics_file           = rna_seq_metrics_file,

            workflow_dot_file              = workflow_dot_file,
    }
}
