version 1.0

######################################################################################
## A workflow that performs single sample genome assembly on PacBio HiFi reads from
## one or more flow cells. The workflow merges multiple samples into a single BAM
## prior to genome assembly and variant calling.
######################################################################################

import "tasks/Utils.wdl" as Utils
import "tasks/Hifiasm.wdl" as HA
import "tasks/CallAssemblyVariants.wdl" as AV
import "tasks/Quast.wdl" as Quast
import "tasks/Finalize.wdl" as FF

workflow PBAssembleWithHifiasm {
    input {
        Array[File] ccs_bams
        Array[File] ccs_pbis

        File ref_map_file
        String participant_name
        String prefix

        String gcs_out_root_dir
    }

    parameter_meta {
        ccs_bams:           "GCS path to unaligned CCS BAM files"
        ccs_pbis:           "GCS path to unaligned CCS BAM file indices"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"
        participant_name:   "name of the participant from whom these samples were obtained"
        prefix:             "prefix for output files"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBAssembleWithHifiasm/~{prefix}"

    call Utils.ComputeGenomeLength { input: fasta = ref_map['fasta'] }

    scatter (ccs_bam in ccs_bams) {
        call Utils.BamToFastq { input: bam = ccs_bam, prefix = basename(ccs_bam, ".bam") }
    }

    # gather across (potential multiple) input CCS BAMs
    if (length(ccs_bams) > 1) {
        call Utils.MergeFastqs as MergeAllFastqs { input: fastqs = BamToFastq.reads_fq }
    }

    File ccs_fq  = select_first([ MergeAllFastqs.merged_fastq, BamToFastq.reads_fq ])

    call HA.Hifiasm {
        input:
            reads = ccs_fq,
            prefix = prefix
    }

    call Quast.Quast {
        input:
            ref = ref_map['fasta'],
            assemblies = [ Hifiasm.fa ]
    }

    call AV.CallAssemblyVariants {
        input:
            asm_fasta = Hifiasm.fa,
            ref_fasta = ref_map['fasta'],
            participant_name = participant_name,
            prefix = prefix + ".hifiasm"
    }

    call FF.FinalizeToFile as FinalizeHifiasmGfa {
        input:
            file    = Hifiasm.gfa,
            outfile = outdir + "/assembly/" + basename(Hifiasm.gfa)
    }

    call FF.FinalizeToFile as FinalizeHifiasmFa {
        input:
            file    = Hifiasm.fa,
            outfile = outdir + "/assembly/" + basename(Hifiasm.fa)
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
        File hifiasm_gfa = FinalizeHifiasmGfa.gcs_path
        File hifiasm_fa = FinalizeHifiasmFa.gcs_path

        File paf = CallAssemblyVariants.paf
        File paftools_vcf = CallAssemblyVariants.paftools_vcf

        File quast_report_html = FinalizeQuastReportHtml.gcs_path
        File quast_report_txt = FinalizeQuastReportTxt.gcs_path

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