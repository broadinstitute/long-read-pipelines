version 1.0

######################################################################################
## A workflow that performs single sample genome assembly on ONT reads from one or
## more flow cells. The workflow merges multiple samples into a single BAM prior to
## genome assembly and variant calling.
######################################################################################

import "tasks/Utils.wdl" as Utils
import "tasks/Flye.wdl" as Flye
import "tasks/Medaka.wdl" as Medaka
import "tasks/CallAssemblyVariants.wdl" as AV
import "tasks/Quast.wdl" as Quast
import "tasks/Finalize.wdl" as FF

workflow ONTAssembleWithFlye {
    input {
        String gcs_fastq_dir

        File ref_map_file

        String medaka_model = "r941_prom_high_g360"

        String participant_name
        String prefix

        String gcs_out_root_dir
    }

    parameter_meta {
        gcs_fastq_dir:       "GCS path to unaligned CCS BAM files"

        ref_map_file:        "table indicating reference sequence and auxillary file locations"

        medaka_model:        "Medaka polishing model name"

        participant_name:    "name of the participant from whom these samples were obtained"
        prefix:              "prefix for output files"

        gcs_out_root_dir:    "GCS bucket to store the reads, variants, and metrics files"
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

    call FF.FinalizeToFile as FinalizeAsmUnpolished {
        input:
            file    = Flye.fa,
            outfile = outdir + "/assembly/" + basename(Flye.fa)
    }

    call FF.FinalizeToFile as FinalizeAsmPolished {
        input:
            file    = MedakaPolish.polished_assembly,
            outfile = outdir + "/assembly/" + basename(MedakaPolish.polished_assembly)
    }

    call FF.FinalizeToFile as FinalizeQuastResults {
        input:
            file    = Quast.results,
            outfile = outdir + "/assembly/" + basename(Quast.results)
    }

    call FF.FinalizeToFile as FinalizeQuastReportHtml {
        input:
            file    = Quast.report_html,
            outfile = outdir + "/assembly/" + basename(Quast.report_html)
    }

    call FF.FinalizeToFile as FinalizeQuastReportTxt{
        input:
            file    = Quast.report_txt,
            outfile = outdir + "/assembly/" + basename(Quast.report_txt)
    }

    output {
        File asm_unpolished = FinalizeAsmUnpolished.gcs_path
        File asm_polished = FinalizeAsmPolished.gcs_path

        File paf = CallAssemblyVariants.paf
        File paftools_vcf = CallAssemblyVariants.paftools_vcf

        File quast_report_html = FinalizeQuastResults.gcs_path
        File quast_report_txt = FinalizeQuastResults.gcs_path

        Int num_contigs = Quast.metrics['#_contigs']
        Int largest_contigs = Quast.metrics['Largest_contig']
        Int total_length = Quast.metrics['Total_length']
        Float genome_fraction_pct = Quast.metrics['Genome_fraction_(%)']
        Float gc_pct = Quast.metrics['GC_(%)']
        Int n50 = Quast.metrics['N50']
        Int ng50 = Quast.metrics['NG50']
        Int nga50 = Quast.metrics['NGA50']
        Int total_aligned_length = Quast.metrics['Total_aligned_length']
        Int largest_alignment = Quast.metrics['Largest_alignment']
        Int unaligned_length = Quast.metrics['Unaligned_length']
        Float duplication_ratio = Quast.metrics['Duplication_ratio']
        Int num_misassemblies = Quast.metrics['#_misassemblies']
        Float num_mismatches_per_100_kbp = Quast.metrics['#_mismatches_per_100_kbp']
        Float num_indels_per_100_kbp = Quast.metrics['#_indels_per_100_kbp']
    }
}