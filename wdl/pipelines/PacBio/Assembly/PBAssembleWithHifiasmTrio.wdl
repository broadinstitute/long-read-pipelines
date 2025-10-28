version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Assembly/HifiasmTrio.wdl" as HAT
import "../../../tasks/QC/Quast.wdl" as QuastEval
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow PBAssembleWithHifiasmTrio {

    meta {
        description: "A workflow that performs single sample genome assembly on PacBio HiFi reads from one or more SMRT cells. The multiple SMRT cells data are merged prior to assembly."
    }
    parameter_meta {
        ccs_fqs:            "GCS path to CCS fastq files"

        participant_name:   "name of the participant from whom these samples were obtained"
        prefix:             "prefix for output files"

        ref_fasta_for_eval: "Reference Fasta used for evaluating "
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        Array[File] ccs_fqs

        String participant_name
        String prefix

        File maternal_fastq_1
        File paternal_fastq_1

        File? ref_fasta_for_eval
        String? quast_extra_args
        File? ont_ultralong_reads
        File? maternal_fastq_2
        File? paternal_fastq_2
        String? telomere_5_prime_sequence

        Boolean haploid = false

        String gcs_out_root_dir
    }

    #########################################################################################
    if (length(ccs_fqs) > 1) {
        call Utils.MergeFastqs as t001_MergeAllFastqs { input: fastqs = ccs_fqs }
    }
    File ccs_fq  = select_first([ t001_MergeAllFastqs.merged_fastq, ccs_fqs[0] ])

    call HAT.HifiasmTrio as t002_HifiasmTrio {
        input:
            reads = ccs_fq,
            ont_ultralong_reads = ont_ultralong_reads,
            maternal_fastq_1 = maternal_fastq_1,
            maternal_fastq_2 = maternal_fastq_2,
            paternal_fastq_1 = paternal_fastq_1,
            paternal_fastq_2 = paternal_fastq_2,
            telomere_5_prime_sequence = telomere_5_prime_sequence,
            haploid = haploid,
            prefix = prefix
    }

    # todo: assumes ploidy 2
    call QuastEval.Quast as t003_primary_h0_h1_quast {
        input:
            ref = ref_fasta_for_eval,
            is_large = true,
            assemblies = [t002_HifiasmTrio.merged_unitigs_fa,
                          t002_HifiasmTrio.maternal_phased_fa,
                          t002_HifiasmTrio.paternal_phased_fa],
            extra_args = quast_extra_args
    }

    call QuastEval.SummarizeQuastReport as t004_primary_h0_h1_quast_summary {
        input: quast_report_txt = t003_primary_h0_h1_quast.report_txt
    }

    #########################################################################################
    # Finalize data
    String workflow_name = "PBAssembleWithHifiasmTrio"

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name + "/~{prefix}"
    String dir = outdir + "/assembly"

    # merged FASTQ
    String dummy = basename(ccs_fq)
    String dummy_b = sub(dummy, ".gz$", "")
    if (dummy != dummy_b) {
        call FF.FinalizeToFile as t004_FinalizeMergedFQ { input: outdir = dir, file = ccs_fq, name = prefix + ".fq.gz" }
    }
    if (dummy == dummy_b) {
        call FF.CompressAndFinalize as t005_CompressAndFinalizeMergedFQ { input: outdir = dir, file = ccs_fq, name = prefix + ".fq.gz" }
    }
    String finalized_merged_fq_path = select_first([t004_FinalizeMergedFQ.gcs_path, t005_CompressAndFinalizeMergedFQ.gcs_path])


    # assembly results themselves
    call FF.CompressAndFinalize as t006_FinalizeHifiasmPrimaryGFA   { input: outdir = dir, file = t002_HifiasmTrio.maternal_phased_gfa }
    call FF.CompressAndFinalize as t007_FinalizeHifiasmPrimaryFA    { input: outdir = dir, file = t002_HifiasmTrio.maternal_phased_fa }

    call FF.CompressAndFinalize as t008_FinalizeHifiasmAlternateGFA   { input: outdir = dir, file = t002_HifiasmTrio.paternal_phased_gfa }
    call FF.CompressAndFinalize as t009_FinalizeHifiasmAlternateFA    { input: outdir = dir, file = t002_HifiasmTrio.paternal_phased_fa }

    call FF.CompressAndFinalize as t010_FinalizeMergedUnitigsGFA   { input: outdir = dir, file = t002_HifiasmTrio.merged_unitigs_gfa }
    call FF.CompressAndFinalize as t011_FinalizeMergedUnitigsFA    { input: outdir = dir, file = t002_HifiasmTrio.merged_unitigs_fa }

    call FF.CompressAndFinalize as t012_FinalizeAllOutputsGZ    { input: outdir = dir, file = t002_HifiasmTrio.all_outputs_gz }

    call FF.FinalizeToFile as t013_FinalizeQuastReportHtml {
        input: outdir = dir, file = t003_primary_h0_h1_quast.report_html
    }
    call FF.FinalizeAndCompress as t014_FinalizeQuastReports {
        input: outdir = dir, files = t003_primary_h0_h1_quast.report_in_various_formats, prefix = prefix + ".quast_reports"
    }
    call FF.FinalizeToFile as t015_FinalizeQuastSummaryAll {
        input: outdir = dir, file = t004_primary_h0_h1_quast_summary.quast_metrics_together
    }
    scatter (report in t004_primary_h0_h1_quast_summary.quast_metrics ) {
        call FF.FinalizeToFile as t016_FinalizeQuastIndividualSummary  { input: outdir = dir, file = report }
    }

    output {
        File merged_fq = finalized_merged_fq_path

        File maternal_phased_gfa  = t006_FinalizeHifiasmPrimaryGFA.gcs_path
        File maternal_phased_fa = t007_FinalizeHifiasmPrimaryFA.gcs_path

        File paternal_phased_gfa  = t008_FinalizeHifiasmAlternateGFA.gcs_path
        File paternal_phased_fa = t009_FinalizeHifiasmAlternateFA.gcs_path

        File merged_unitigs_gfa = t010_FinalizeMergedUnitigsGFA.gcs_path
        File merged_unitigs_fa = t011_FinalizeMergedUnitigsFA.gcs_path

        File all_outputs_gz = t012_FinalizeAllOutputsGZ.gcs_path

        File? quast_report_html = t013_FinalizeQuastReportHtml.gcs_path
        File? quast_report_in_various_formats = t014_FinalizeQuastReports.gcs_path

        File? quast_summary_on_all = t015_FinalizeQuastSummaryAll.gcs_path

        File? quast_summary_on_merged_unitigs = t016_FinalizeQuastIndividualSummary.gcs_path[0]
        File? quast_summary_on_maternal_haplotype = t016_FinalizeQuastIndividualSummary.gcs_path[1]
        File? quast_summary_on_paternal_haplotype = t016_FinalizeQuastIndividualSummary.gcs_path[2]
    }
}
