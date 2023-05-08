version 1.0

import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/Assembly/Hifiasm.wdl" as HA
import "../../../tasks/QC/Quast.wdl" as QuastEval

workflow PBAssembleWithHifiasm {

    meta {
        description: "A workflow that performs single, diploid sample genome assembly on PacBio HiFi reads from one or more SMRT cells. The multiple SMRT cells data are merged prior to assembly."
    }
    parameter_meta {
        ccs_fqs:            "GCS path to CCS fastq files"
        prefix:             "prefix for output files"

        ref_fasta_for_eval: "Reference Fasta used for evaluating "
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"

        hifiasm_primary_gfa: "primary assembly GFA output (gzipped)"

        hifiasm_alternate_tigs: "alternative assembly FASTA output (block gzipped)"

        hifiasm_haploGFAs:   "path to folder hosting GFA files (gzipped) for the haplotype-resolved assemblies"
        hifiasm_haplotigs:   "path to folder hosting FASTA files (block gzipped) for the haplotype-resolved assemblies"

        quast_report_html:    "QUAST report on [primary, H0, H1] assemblies FASTA, in HTML format"
        quast_summary_on_all: "QUAST summary on [primary, H0, H1] assemblies FASTA"
    }

    input {
        Array[File] ccs_fqs

        String prefix

        File? ref_fasta_for_eval

        String gcs_out_root_dir

        Array[String] gcp_zones = ['us-central1-a', 'us-central1-b', 'us-central1-c', 'us-central1-f']
    }

    #########################################################################################
    if (length(ccs_fqs) > 1) {
        call Utils.MergeFastqs as MergeAllFastqs { input: fastqs = ccs_fqs }
    }
    File ccs_fq  = select_first([ MergeAllFastqs.merged_fastq, ccs_fqs[0] ])

    #########################################################################################
    call GU.CollapseArrayOfStrings as get_zones {input: input_array = gcp_zones, joiner = " "}
    String wdl_parsable_zones = get_zones.collapsed

    call HA.Hifiasm {
        input:
            reads = ccs_fq,
            prefix = prefix,
            zones = wdl_parsable_zones
    }

    # todo: assumes ploidy 2
    call QuastEval.Quast as primary_h0_h1_quast {
        input:
            ref = ref_fasta_for_eval,
            is_large = true,
            assemblies = [Hifiasm.primary_tigs,
                          Hifiasm.hap1_tigs,
                          Hifiasm.hap2_tigs]
    }

    call QuastEval.SummarizeQuastReport as primary_h0_h1_quast_summary {
        input: quast_report_txt = primary_h0_h1_quast.report_txt
    }

    #########################################################################################
    # Finalize data
    String workflow_name = "PBAssembleWithHifiasm"

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name + "/~{prefix}"

    ##########
    # merged FASTQ
    String dummy = basename(ccs_fq)
    String dummy_b = sub(dummy, ".gz$", "")
    if (dummy != dummy_b) {
        call FF.FinalizeToFile as FinalizeMergedFQ { input: outdir = outdir, file = ccs_fq, name = prefix + ".fq.gz" }
    }
    if (dummy == dummy_b) {
        call FF.CompressAndFinalize as CompressAndFinalizeMergedFQ { input: outdir = outdir, file = ccs_fq, name = prefix + ".fq.gz" }
    }
    String finalized_merged_fq_path = select_first([FinalizeMergedFQ.gcs_path, CompressAndFinalizeMergedFQ.gcs_path])

    ##########
    # assembly results themselves
    String asm_dir = outdir + "/assembly"

    # primary/alt
    call FF.CompressAndFinalize as FinalizeHifiasmPrimaryGFA   { input: outdir = asm_dir, file = Hifiasm.primary_gfa }
    call FF.FinalizeToFile      as FinalizeHifiasmPrimaryFA    { input: outdir = asm_dir, file = Hifiasm.primary_tigs }
    call FF.FinalizeToFile      as FinalizeHifiasmPrimaryGzi   { input: outdir = asm_dir, file = Hifiasm.primary_tigs_gzi }

    call FF.CompressAndFinalize as FinalizeHifiasmAlternateGFA   { input: outdir = asm_dir, file = Hifiasm.alternate_gfa }
    call FF.FinalizeToFile      as FinalizeHifiasmAlternateFA    { input: outdir = asm_dir, file = Hifiasm.alternate_tigs }
    call FF.FinalizeToFile      as FinalizeHifiasmAlternateGzi   { input: outdir = asm_dir, file = Hifiasm.alternate_tigs_gzi }

    # H1/H2
    call FF.FinalizeToFile as FinalizeHifiasmHapOneGFA   { input: outdir = asm_dir, file = Hifiasm.hap1_gfa }
    call FF.FinalizeToFile as FinalizeHifiasmHapOneFASTA { input: outdir = asm_dir, file = Hifiasm.hap1_tigs }
    call FF.FinalizeToFile as FinalizeHifiasmHapOneFaGzi { input: outdir = asm_dir, file = Hifiasm.hap1_tig_gzi }
    call FF.FinalizeToFile as FinalizeHifiasmHapTwoGFA   { input: outdir = asm_dir, file = Hifiasm.hap2_gfa }
    call FF.FinalizeToFile as FinalizeHifiasmHapTwoFASTA { input: outdir = asm_dir, file = Hifiasm.hap2_tigs }
    call FF.FinalizeToFile as FinalizeHifiasmHapTwoFaGzi { input: outdir = asm_dir, file = Hifiasm.hap2_tig_gzi }

    call FF.FinalizeToFile as FinalizeHifiasmLog { input:
        outdir = asm_dir, file = Hifiasm.log_in_hap_mode, name = "~{prefix}.hifiasm-hapmode.log"
    }
    call FF.FinalizeToFile as FinalizeHifiasmResourceUsagesVisual { input:
        outdir = asm_dir, file = Hifiasm.resource_usage_visual_in_hap_mode
    }

    ##########
    # quast stuff
    String quast_dir = outdir + "/quast"
    call FF.FinalizeToFile as FinalizeQuastReportHtml {
        input: outdir = quast_dir, file = primary_h0_h1_quast.report_html, name = '~{prefix}.quast-report.html'
    }
    call FF.FinalizeToFile as FinalizeQuastSummaryAll {
        input: outdir = quast_dir, file = primary_h0_h1_quast_summary.quast_metrics_together, name = '~{prefix}.quast-summary.txt'
    }
    call FF.FinalizeToDir as FinalizeQuastReports { input:
        outdir = quast_dir + "/misc/",

        files = flatten([primary_h0_h1_quast.plots,
                         primary_h0_h1_quast.report_in_various_formats,
                         primary_h0_h1_quast_summary.quast_metrics])
    }
    if (defined(primary_h0_h1_quast.contigs_reports)) {
        call FF.FinalizeToFile as FinalizeQuastContigsReport { input: outdir = quast_dir, file = select_first([primary_h0_h1_quast.contigs_reports])}
    }

    ###########################################################
    call GU.GetTodayDate as today {}

    ###########################################################
    output {
        String last_processing_date = today.yyyy_mm_dd

        ########################################
        File merged_fq = finalized_merged_fq_path

        ########################################
        Map[String, File] hifiasm_hap_outputs = {
            "HapOne_GFA": FinalizeHifiasmHapOneGFA.gcs_path,
            "HapTwo_GFA": FinalizeHifiasmHapTwoGFA.gcs_path,
            "HapOne_FASTA": FinalizeHifiasmHapOneFASTA.gcs_path,
            "HapTwo_FASTA": FinalizeHifiasmHapTwoFASTA.gcs_path,
            "HapOne_FA_GZI": FinalizeHifiasmHapOneFaGzi.gcs_path,
            "HapTWO_FA_GZI": FinalizeHifiasmHapTwoFaGzi.gcs_path,
            "resource_use_visual": FinalizeHifiasmResourceUsagesVisual.gcs_path,
            "runtime_log": FinalizeHifiasmLog.gcs_path,
        }

        Map[String, File] hifiasm_primary_alt_outputs = {
            "primary_gfa": FinalizeHifiasmPrimaryGFA.gcs_path,
            "primary_fasta": FinalizeHifiasmPrimaryFA.gcs_path,
            "primary_fa_gzi": FinalizeHifiasmPrimaryGzi.gcs_path,
            "alternate_gfa": FinalizeHifiasmAlternateGFA.gcs_path,
            "alternate_fasta": FinalizeHifiasmAlternateFA.gcs_path,
            "alternate_fa_gzi": FinalizeHifiasmAlternateGzi.gcs_path,
        }

        ########################################
        File quast_report_html = FinalizeQuastReportHtml.gcs_path
        File quast_summary_on_all = FinalizeQuastSummaryAll.gcs_path

        Map[String, String] quast_2ndary_outputs = {
            "quast_2ndary_outputs": FinalizeQuastReports.gcs_dir,
            "quast_contigs_report": if (defined(primary_h0_h1_quast.contigs_reports)) then select_first([FinalizeQuastContigsReport.gcs_path]) else "None"
        }
    }
}
