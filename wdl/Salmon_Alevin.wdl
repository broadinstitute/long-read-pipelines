version 1.0

import "tasks/TranscriptAnalysis/Salmon_Tasks.wdl" as SALMON

workflow Salmon_Alevin {

    meta {
        description : "Quantify transcripts from single cell RNA transcript reads using SALMON Alevin."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File mates1_fastq
        File mates2_fastq
        File cb_whitelist

        File transcript_fasta = "gs://broad-dsde-methods-long-reads/resources/gencode_v34/gencode.v34.pc_transcripts.fa"
        File tgmap = "gs://broad-dsde-methods-long-reads/resources/gencode_v34/gencode.v34.pc_transcripts.tgmap.tsv"
        Int index_kmer_size = 31
        Boolean is_gencode_transcriptome = true

        String library_type = "ISR"
        Int end = 5
        Int umi_length = 10
        Int barcode_length = 16

        String prefix = "alevin_output"
    }

    parameter_meta {
        mates1_fastq : "FASTQ file containing the mates1 information for Alevin (CBC + UMI for each read)."
        mates2_fastq : "FASTQ file containing the mates2 information for Alevin (transcript sequence for each read)."
        cb_whitelist : "Cell barcode whitelist file for Alevin.  One cell barcode per line."

        transcript_fasta : "[optional] FASTA file containing all transcript sequences to quantify.  (Default: gs://broad-dsde-methods-long-reads/resources/gencode_v34/gencode.v34.pc_transcripts.fa)"
        tgmap : "Transcript / Gene map for Alevin.  Unheadered TSV with two columns - Transcript | Gene.  (Default: gs://broad-dsde-methods-long-reads/resources/gencode_v34/gencode.v34.pc_transcripts.tgmap.tsv)"
        index_kmer_size : "[optional] K-mer size to use for the SALMON index of the given transcript_fasta.  (Default: 31)"
        is_gencode_transcriptome : "[optional] True if the given transcript_fasta is a gencode file.  False otherwise. (Default: True)"

        library_type : "[optional] SALMON library type (Default: ISR)."
        end : "[optional] End of the library prep.  (Default: 5)"
        umi_length : "[optional] Length of the unique molecular identifier sequence in the library in bases.  (Default: 10)"
        barcode_length : "[optional] Length of the unique molecular identifier sequence in the library in bases.  (Default: 16)"

        prefix : "[optional] Prefix to use for the created output files.  (Default: alevin_output)"
    }

    call SALMON.Alevin as alevin {
        input:
            mates1_fastq = mates1_fastq,
            mates2_fastq = mates2_fastq,
            cb_whitelist = cb_whitelist,
            transcript_fasta = transcript_fasta,
            tgmap = tgmap,
            index_kmer_size = index_kmer_size,
            is_gencode_transcriptome = is_gencode_transcriptome,
            library_type = library_type,
            end = end,
            umi_length = umi_length,
            barcode_length = barcode_length,
            prefix = prefix,
    }

    # ------------------------------------------------
    # Outputs:
    output {
        File output_tar_gz = alevin.output_tar_gz
        File output_fastq = alevin.output_fastq
    }
}