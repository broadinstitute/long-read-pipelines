version 1.0

import "Flair_Tasks.wdl" as FLAIR

workflow RunFlairQuant {

    meta {
        description : "Quantify transcript isoforms using FLAIR quantify."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        # ------------------------------------------------
        # Input args:
        # Required:

        File fasta_tar_gz
        File transcript_isoforms_fasta
        String out_base_name = "reads"

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    parameter_meta {
        fasta_tar_gz : "tar.gz file contianing FASTA files, each with reads from the same sample / cell barcode to quantify for transcripts.  Corresponds to `fasta_tar_gz_out` from Preprocessing_Tasks.SplitBamBySampleAndCellBarcodeTask."
        transcript_isoforms_fasta : "FASTA file containing isoforms sequences to quantify."
        out_base_name : "[optional] Base name for the resulting output file.  Ideally should be a unique identifier for this dataset."

        mem_gb : "[optional] Amount of memory to give to the machine running each task in this workflow."
        preemptible_attempts : "[optional] Number of times to allow each task in this workflow to be preempted."
        disk_space_gb : "[optional] Amount of storage disk space (in Gb) to give to each machine running each task in this workflow."
        cpu : "[optional] Number of CPU cores to give to each machine running each task in this workflow."
        boot_disk_size_gb : "[optional] Amount of boot disk space (in Gb) to give to each machine running each task in this workflow."
    }

    call FLAIR.FlairQuant as quant {
        input:
            fasta_tar_gz              = fasta_tar_gz,
            transcript_isoforms_fasta = transcript_isoforms_fasta
    }

    # ------------------------------------------------
    # Outputs:
    output {
      # Default output file name:
      File count_matrix         = quant.count_matrix

      File timing_info          = quant.timing_info
      File memory_log           = quant.memory_log
    }
}
