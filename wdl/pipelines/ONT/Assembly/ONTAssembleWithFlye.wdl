version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Assembly/Flye.wdl" as Flye
import "../../../tasks/Preprocessing/Medaka.wdl" as Medaka
import "../../../tasks/VariantCalling/CallAssemblyVariants.wdl" as AV
import "../../../tasks/QC/Quast.wdl" as Quast
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ONTAssembleWithFlye {
    meta {
        description: "Perform single sample genome assembly on ONT reads from one or more flow cells. The workflow merges multiple samples into a single BAM prior to genome assembly and variant calling."
    }
    parameter_meta {
        gcs_fastq_dir:       "GCS path to unaligned CCS BAM files"

        ref_map_file:        "table indicating reference sequence and auxillary file locations"

        medaka_model:        "Medaka polishing model name"

        participant_name:    "name of the participant from whom these samples were obtained"
        prefix:              "prefix for output files"

        gcs_out_root_dir:    "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        String gcs_fastq_dir

        File ref_map_file

        String medaka_model = "r941_prom_high_g360"

        String participant_name
        String prefix

        String gcs_out_root_dir
    }


    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTAssembleWithFlye/~{prefix}"

    call Utils.ComputeGenomeLength { input: fasta = ref_map['fasta'] }

    call Utils.ListFilesOfType { input: gcs_dir = gcs_fastq_dir, suffixes = [".fastq", ".fq", ".fastq.gz", ".fq.gz"] }
    call Utils.MergeFastqs { input: fastqs = ListFilesOfType.files }

    call Flye.Flye {
        input:
            reads = MergeFastqs.merged_fastq,
            genome_size = ComputeGenomeLength.length,
            prefix = prefix,
    }

    call Medaka.MedakaPolish {
        input:
            basecalled_reads = MergeFastqs.merged_fastq,
            draft_assembly = Flye.fa,
            model = medaka_model,
            prefix = basename(Flye.fa, ".fasta") + ".consensus",
            n_rounds = 3
    }

    call Quast.Quast {
        input:
            ref = ref_map['fasta'],
            assemblies = [ MedakaPolish.polished_assembly ]
    }

    call AV.CallAssemblyVariants {
        input:
            asm_fasta = MedakaPolish.polished_assembly,
            ref_fasta = ref_map['fasta'],
            participant_name = participant_name,
            prefix = prefix + ".canu"
    }

    # Finalize data
    String dir = outdir + "/assembly"

    call FF.FinalizeToFile as FinalizeAsmUnpolished { input: outdir = dir, file = Flye.fa }
    call FF.FinalizeToFile as FinalizeAsmPolished { input: outdir = dir, file = MedakaPolish.polished_assembly }
    call FF.FinalizeToFile as FinalizeQuastReportHtml { input: outdir = dir, file = Quast.report_html }
    call FF.FinalizeToFile as FinalizeQuastReportTxt { input: outdir = dir, file = Quast.report_txt }
    call FF.FinalizeToFile as FinalizePaf { input: outdir = dir, file = CallAssemblyVariants.paf }
    call FF.FinalizeToFile as FinalizePafToolsVcf { input: outdir = dir, file = CallAssemblyVariants.paftools_vcf }

    call Quast.SummarizeQuastReport as summaryQ {input: quast_report_txt = Quast.report_txt}
    Map[String, String] q_metrics = read_map(summaryQ.quast_metrics[0])

    output {
        File asm_unpolished = FinalizeAsmUnpolished.gcs_path
        File asm_polished = FinalizeAsmPolished.gcs_path

        File paf = FinalizePaf.gcs_path
        File paftools_vcf = FinalizePafToolsVcf.gcs_path

        File quast_report_html = FinalizeQuastReportHtml.gcs_path
        File quast_report_txt = FinalizeQuastReportTxt.gcs_path

        Map[String, String] quast_summary = q_metrics
    }
}