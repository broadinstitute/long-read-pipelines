version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/StringTie2.wdl"
import "tasks/Finalize.wdl" as FF

workflow AnnotateTranscriptomeWithGuide {
    input {
        File aligned_bam
        File aligned_bai
        File gtf
        File ref_map_file

        String participant_name

        String gcs_out_root_dir
    }

    parameter_meta {
        aligned_bam:        "GCS path to aligned BAM file"
        aligned_bai:        "GCS path to aligned BAM file index"

        gtf:                "guide GTF (e.g. Gencode38)"
        ref_map_file:       "table indicating reference sequence and auxillary file locations"
        participant_name:   "name of the participant from whom these samples were obtained"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/AnnotateTranscriptomeWithGuide/~{participant_name}"

    call StringTie2.Quantify as QuantifyWithoutRetainedIntrons {
        input:
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai,
            gtf = gtf,
            keep_retained_introns = false,
            prefix = participant_name + ".without_retained_introns"
    }

    call StringTie2.Quantify as QuantifyWithRetainedIntrons {
        input:
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai,
            gtf = gtf,
            keep_retained_introns = true,
            prefix = participant_name + ".with_retained_introns"
    }

    call StringTie2.ExtractTranscriptSequences as ExtractTranscriptsWithoutRetainedIntrons {
        input:
            ref_fasta = ref_map['fasta'],
            ref_fasta_fai = ref_map['fai'],
            gtf = QuantifyWithoutRetainedIntrons.st_gtf
    }

    call StringTie2.ExtractTranscriptSequences as ExtractTranscriptsWithRetainedIntrons {
        input:
            ref_fasta = ref_map['fasta'],
            ref_fasta_fai = ref_map['fai'],
            gtf = QuantifyWithRetainedIntrons.st_gtf
    }

    # Finalize
    String idir = outdir + "/transcriptome/with_retained_introns"

    call FF.FinalizeToFile as FinalizeGTFWithIntrons {
        input:
            outdir = idir,
            file = QuantifyWithRetainedIntrons.st_gtf,
            name = "~{participant_name}.with_retained_introns.gtf"
    }

    call FF.FinalizeToFile as FinalizeFaWithIntrons {
        input:
            outdir = idir,
            file = ExtractTranscriptsWithRetainedIntrons.transcripts_fa,
            name = "~{participant_name}.with_retained_introns.fa"
    }

    String odir = outdir + "/transcriptome/without_retained_introns"

    call FF.FinalizeToFile as FinalizeGTFWithoutIntrons {
        input:
            outdir = odir,
            file = QuantifyWithoutRetainedIntrons.st_gtf,
            name = "~{participant_name}.without_retained_introns.gtf"
    }

    call FF.FinalizeToFile as FinalizeFaWithoutIntrons {
        input:
            outdir = odir,
            file = ExtractTranscriptsWithoutRetainedIntrons.transcripts_fa,
            name = "~{participant_name}.without_retained_introns.fa"
    }

    output {
        File with_retained_introns_gtf = FinalizeGTFWithIntrons.gcs_path
        File with_retained_introns_fa = FinalizeFaWithIntrons.gcs_path

        File without_retained_introns_gtf = FinalizeGTFWithoutIntrons.gcs_path
        File without_retained_introns_fa = FinalizeFaWithoutIntrons.gcs_path
    }
}
