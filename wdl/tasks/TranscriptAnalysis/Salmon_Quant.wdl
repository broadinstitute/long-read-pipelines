version 1.0

import "Salmon_Tasks.wdl" as SALMON

workflow RunSalmonOnCellBarcodedReads {

    meta {
        description : "Quantify transcripts from RNA transcript reads using SALMON."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        Array[File] cell_fastas
        File salmon_index_tar_gz = "gs://broad-dsde-methods-long-reads/resources/gencode_v34/gencode.v34.pc_transcripts_index_k31.tar.gz"
    }

    parameter_meta {
        cell_fastas : "Array of FASTA/FASTQ files."
        salmon_index_tar_gz : "[optional] SALMON index file corresponding to the transcripts FASTA file used to quantify the transcripts in the reads (TAR.GZ format)."
    }

    scatter (reads_fasta in cell_fastas) {

        call SALMON.RunSalmonQuantTask as salmon_quant {
            input:
                reads_fasta = reads_fasta,
                salmon_index_tar_gz = salmon_index_tar_gz,
                extra_args = "--fldMean 751 --fldSD 460 --minAssignedFrags 1"

        }
    }

    call SALMON.ConvertQuantFilesToCountMatrix as make_count_matrix {
        input:
            quant_files = salmon_quant.quant_file
    }

    # ------------------------------------------------
    # Outputs:
    output {
      Array[File] quant_files       = salmon_quant.quant_file
      Array[File] cmd_infos         = salmon_quant.cmd_info
      Array[File] lib_format_counts = salmon_quant.lib_format_counts
      Array[File] ambig_infos       = salmon_quant.ambig_info
      Array[File] eq_classes        = salmon_quant.eq_classes
      Array[File] expected_biases   = salmon_quant.expected_bias
      Array[File] flds              = salmon_quant.fld
      Array[File] meta_infos        = salmon_quant.meta_info
      Array[File] observed_biases   = salmon_quant.observed_bias
      Array[File] observed_bias_3ps = salmon_quant.observed_bias_3p
      Array[File] logs              = salmon_quant.log

      File salmon_count_matrix_tsv  = make_count_matrix.count_matrix_tsv
      File salmon_count_matrix_h5ad = make_count_matrix.count_matrix_h5ad
    }
}