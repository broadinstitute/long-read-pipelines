version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/AsmUtils.wdl" as AU
import "tasks/Canu.wdl" as Canu
import "tasks/Medaka.wdl" as Medaka
import "tasks/Quast.wdl" as Quast
import "tasks/Finalize.wdl" as FF

workflow ONTAssembleWithCanu {
    input {
        String gcs_fastq_dir

        File ref_map_file

        Float correct_error_rate = 0.15
        Float trim_error_rate = 0.15
        Float assemble_error_rate = 0.15
        String medaka_model = "r941_prom_high_g360"

        String prefix

        String gcs_out_root_dir
    }

    Map[String, String] ref_map = read_map(ref_map_file)
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTAssembleWithCanu/~{prefix}"

    call Utils.ComputeGenomeLength { input: fasta = ref_map['fasta'] }

    call Utils.ListFilesOfType { input: gcs_dir = gcs_fastq_dir, suffixes = [".fastq", ".fq", ".fastq.gz", ".fq.gz"] }
    call AU.MergeFastqs { input: fastqs = ListFilesOfType.files }

    call Canu.Canu {
        input:
            reads = MergeFastqs.merged_fastq,
            prefix = prefix,
            genome_size = ceil(ComputeGenomeLength.length/1000000.0),
            correct_error_rate = correct_error_rate,
            trim_error_rate = trim_error_rate,
            assemble_error_rate = assemble_error_rate
    }

    call Medaka.MedakaPolish {
        input:
            basecalled_reads = MergeFastqs.merged_fastq,
            draft_assembly = Canu.fa,
            model = medaka_model,
            prefix = basename(Canu.fa, ".fasta") + ".consensus",
            n_rounds = 3
    }

    call Quast.Quast {
        input:
            ref = ref_map['fasta'],
            assemblies = [ MedakaPolish.polished_assembly ]
    }

    call FF.FinalizeToFile as FinalizeAsmUnpolished {
        input:
            file    = Canu.fa,
            outfile = outdir + "/assembly/" + basename(Canu.fa)
    }

    call FF.FinalizeToFile as FinalizeAsmPolished {
        input:
            file    = MedakaPolish.polished_assembly,
            outfile = outdir + "/assembly/" + basename(MedakaPolish.polished_assembly)
    }

    call FF.FinalizeToFile as FinalizeQuastReport {
        input:
            file    = Quast.report,
            outfile = outdir + "/assembly/" + basename(Quast.report)
    }

    call FF.FinalizeToFile as FinalizeQuastResults {
        input:
            file    = Quast.results,
            outfile = outdir + "/assembly/" + basename(Quast.results)
    }

    output {
        File asm_unpolished = FinalizeAsmUnpolished.gcs_path
        File asm_polished = FinalizeAsmPolished.gcs_path

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