version 1.0

import "../../../tasks/Utility/StringTie2.wdl"
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow AnnotateTranscriptomeWithGuide {

    meta {
        description: "Annotate a transcriptome with a guide GTF"
    }
    parameter_meta {
        aligned_bam:           "GCS path to aligned BAM file"
        aligned_bai:           "GCS path to aligned BAM file index"
        ref_map_file:          "table indicating reference sequence and auxillary file locations"

        gtf:                   "guide GTF (e.g. Gencode38)"
        participant_name:      "name of the participant from whom these samples were obtained"
        keep_retained_introns: "keep apparently retained introns in new transcriptome annotation"

        gcs_out_root_dir:      "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        File aligned_bam
        File aligned_bai
        File ref_map_file

        File gtf
        String participant_name

        Boolean keep_retained_introns = false

        String gcs_out_root_dir
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String prefix = if keep_retained_introns then "with_retained_introns" else "without_retained_introns"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/AnnotateTranscriptomeWithGuide/~{participant_name}/transcriptome/~{prefix}"

    call StringTie2.Quantify {
        input:
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai,
            gtf = gtf,
            keep_retained_introns = keep_retained_introns,
            prefix = "~{participant_name}.~{prefix}"
    }

    call StringTie2.ExtractTranscriptSequences {
        input:
            ref_fasta = ref_map['fasta'],
            ref_fasta_fai = ref_map['fai'],
            gtf = Quantify.st_gtf,
            prefix = "~{participant_name}.~{prefix}"
    }

    call StringTie2.CompareTranscriptomes {
        input:
            guide_gtf = gtf,
            new_gtf = Quantify.st_gtf,
            prefix = "~{participant_name}.~{prefix}"
    }

    # Finalize
    call FF.FinalizeToFile as FinalizeGTF { input: outdir = outdir, file = Quantify.st_gtf }

    call FF.FinalizeToFile as FinalizeFa { input: outdir = outdir, file = ExtractTranscriptSequences.transcripts_fa }
    call FF.FinalizeToFile as FinalizeFai { input: outdir = outdir, file = ExtractTranscriptSequences.transcripts_fai }
    call FF.FinalizeToFile as FinalizeDict { input: outdir = outdir, file = ExtractTranscriptSequences.transcripts_dict }

    call FF.FinalizeToFile as FinalizeAnnotatedGTF { input: outdir = outdir, file = CompareTranscriptomes.annotated_gtf }
    call FF.FinalizeToFile as FinalizeLoci { input: outdir = outdir, file = CompareTranscriptomes.loci }
    call FF.FinalizeToFile as FinalizeStats { input: outdir = outdir, file = CompareTranscriptomes.stats, name = basename(CompareTranscriptomes.stats) + ".stats.txt" }
    call FF.FinalizeToFile as FinalizeTracking { input: outdir = outdir, file = CompareTranscriptomes.tracking }
    call FF.FinalizeToFile as FinalizeRefMap { input: outdir = outdir, file = CompareTranscriptomes.refmap }
    call FF.FinalizeToFile as FinalizeTMap { input: outdir = outdir, file = CompareTranscriptomes.tmap }

    output {
        File transcriptome_gtf = FinalizeGTF.gcs_path

        File transcriptome_fa = FinalizeFa.gcs_path
        File transcriptome_fai = FinalizeFai.gcs_path
        File transcriptome_dict = FinalizeDict.gcs_path

        File comp_annotated_gtf = FinalizeAnnotatedGTF.gcs_path
        File comp_loci = FinalizeLoci.gcs_path
        File comp_stats = FinalizeStats.gcs_path
        File comp_tracking = FinalizeTracking.gcs_path
        File comp_refmap = FinalizeRefMap.gcs_path
        File comp_tmap = FinalizeTMap.gcs_path
    }
}
