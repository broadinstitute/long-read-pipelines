version 1.0

######################################################################################
## A workflow that performs single sample genome assembly on PacBio HiFi reads from
## one or more flow cells. The workflow merges multiple samples into a single BAM
## prior to genome assembly and variant calling.
######################################################################################

import "tasks/Utils.wdl" as Utils
import "tasks/Hifiasm.wdl" as HA
import "tasks/CallAssemblyVariants.wdl" as AV
import "tasks/Quast.wdl" as QuastEval
import "tasks/Finalize.wdl" as FF

workflow PBAssembleWithHifiasm {
    input {
        Array[File] ccs_fqs

        String participant_name
        String prefix

        File ref_map_file

        String gcs_out_root_dir
    }

    parameter_meta {
        ccs_fqs:            "GCS path to CCS fastq files"

        participant_name:   "name of the participant from whom these samples were obtained"
        prefix:             "prefix for output files"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBAssembleWithHifiasm/~{prefix}"

    if (length(ccs_fqs) > 1) {
        call Utils.MergeFastqs as MergeAllFastqs { input: fastqs = ccs_fqs }
    }

    File ccs_fq  = select_first([ MergeAllFastqs.merged_fastq, ccs_fqs[0] ])

    call HA.Hifiasm {
        input:
            reads = ccs_fq,
            prefix = prefix
    }

    call QuastEval.Quast {
        input:
            assemblies = [Hifiasm.primary_tigs,
                          Hifiasm.phased_tigs[0],
                          Hifiasm.phased_tigs[1]]
    }

    call AV.CallAssemblyVariants {
        input:
            asm_fasta = Hifiasm.primary_tigs,
            ref_fasta = ref_map['fasta'],
            participant_name = participant_name,
            prefix = prefix + ".hifiasm"
    }

    # Finalize data
    String dir = outdir + "/assembly"

    call FF.CompressAndFinalize as FinalizeHifiasmPrimaryGfa   { input: outdir = dir, file = Hifiasm.primary_gfa,  gzip_compress = true }
    call FF.CompressAndFinalize as FinalizeHifiasmPrimaryFa    { input: outdir = dir, file = Hifiasm.primary_tigs, gzip_compress = true }

    call FF.CompressAndFinalize as FinalizeHifiasmAlternateGfa   { input: outdir = dir, file = Hifiasm.alternate_gfa,  gzip_compress = true }
    call FF.CompressAndFinalize as FinalizeHifiasmAlternateFa    { input: outdir = dir, file = Hifiasm.alternate_tigs, gzip_compress = true }

    # we prefer the primary files generated from the primary VS alternate, but if we decide to get that from the haplotig mode, simply finalize them

    call FF.FinalizeAndCompress as FinalizeHifiasmHapGfas  { input: outdir = dir, files = Hifiasm.phased_gfas, prefix = prefix + ".haploGFAs" }
    call FF.FinalizeAndCompress as FinalizeHifiasmHapTigs  { input: outdir = dir, files = Hifiasm.phased_tigs, prefix = prefix + ".haploTigs" }

    call FF.FinalizeToFile as FinalizeQuastReportHtml { input: outdir = dir, file = Quast.report_html }
    call FF.FinalizeToFile as FinalizeQuastReportTxt  { input: outdir = dir, file = Quast.report_txt }
    scatter (report in Quast.quast_metrics) {
        call FF.FinalizeToFile as FinalizeQuastIndividualSummary  { input: outdir = dir, file = report }
    }
    Map[String, String] metrics = read_map(Quast.quast_metrics[0])  # assuming the 0-th is the special

    call FF.FinalizeToFile as FinalizePaf             { input: outdir = dir, file = CallAssemblyVariants.paf }
    call FF.FinalizeToFile as FinalizePafToolsVcf     { input: outdir = dir, file = CallAssemblyVariants.paftools_vcf }

    output {
        File hifiasm_primary_gfa  = FinalizeHifiasmPrimaryGfa.gcs_path
        File hifiasm_primary_tigs = FinalizeHifiasmPrimaryFa.gcs_path

        File hifiasm_haplogfas = FinalizeHifiasmHapGfas.gcs_path
        File hifiasm_haplotigs = FinalizeHifiasmHapTigs.gcs_path

        File hifiasm_alternate_gfa  = FinalizeHifiasmAlternateGfa.gcs_path
        File hifiasm_alternate_tigs = FinalizeHifiasmAlternateFa.gcs_path

        File paf = FinalizePaf.gcs_path
        File paftools_vcf = FinalizePafToolsVcf.gcs_path

        File quast_report_html = FinalizeQuastReportHtml.gcs_path
        File quast_report_txt  = FinalizeQuastReportTxt.gcs_path
        Array[File] individual_quast_summaries = FinalizeQuastIndividualSummary.gcs_path

    #     Int num_contigs = metrics['#_contigs']
    #     Int largest_contigs = metrics['Largest_contig']
    #     String total_length = metrics['Total_length']
    #     Float gc_pct = metrics['GC_(%)']
    #     Int n50 = metrics['N50']
    #     Int n75 = metrics['N75']
    #     Int l50 = metrics['L50']
    #     Int l75 = metrics['L75']

    #    Float genome_fraction_pct = metrics['Genome_fraction_(%)']
    #    Int ng50 = metrics['NG50']
    #    Int nga50 = metrics['NGA50']
    #    Int total_aligned_length = metrics['Total_aligned_length']
    #    Int largest_alignment = metrics['Largest_alignment']
    #    Int unaligned_length = metrics['Unaligned_length']
    #    Float duplication_ratio = metrics['Duplication_ratio']
    #    Int num_misassemblies = metrics['#_misassemblies']
    #    Float num_mismatches_per_100_kbp = metrics['#_mismatches_per_100_kbp']
    #    Float num_indels_per_100_kbp = metrics['#_indels_per_100_kbp']
    }
}