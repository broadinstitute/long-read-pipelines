version 1.0

import "../../../tasks/Utility/PBUtils.wdl" as PB
import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow PBCCSDemultiplex {

    meta {
        description: "A workflow that performs demultiplexing on PacBio CCS data with PacBio's lima tool."
    }
    parameter_meta {
        bam:              "GCS path to CCS BAM file"
        pbi:              "GCS path to CCS BAM .pbi index"
        barcode_file:     "GCS path to the fasta file that specifies the expected set of multiplexing barcodes"

        gcs_out_root_dir: "GCS bucket to store the demultiplexed files and metrics files"
    }

    input {
        File bam
        File pbi
        File barcode_file

        String gcs_out_root_dir
    }

    String basename = basename(bam, ".bam")
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCCSDemultiplex/~{basename}"

    String rdir = outdir + "/reads"
    String mdir = outdir + "/metrics/combined/lima"
    String fdir = outdir + "/figures"

    # demultiplex CCS-ed BAM
    call PB.Demultiplex {
        input:
            bam = bam,
            prefix = basename(bam, ".bam"),
            barcode_file = barcode_file,
            ccs = true,
            guess = 75,
            guess_min_count = 1,
            dump_removed = false,
            split_bam_named = true
    }

    # make reports on demultiplexing
    call PB.MakeSummarizedDemultiplexingReport as SummarizedDemuxReportPNG { input: report = Demultiplex.report }
    call PB.MakeDetailedDemultiplexingReport as DetailedDemuxReportPNG { input: report = Demultiplex.report, type="png" }

    scatter (demux_bam in Demultiplex.demux_bams) {
        call PB.PBIndex as IndexAlignedReads { input: bam = demux_bam }

        # Finalize
        call FF.FinalizeToFile as FinalizeBam { input: outdir = rdir, file = demux_bam }
        call FF.FinalizeToFile as FinalizePbi { input: outdir = rdir, file = IndexAlignedReads.pbi }
    }

    call FF.FinalizeToFile as FinalizeDemultiplexCounts { input: outdir = mdir, file = Demultiplex.counts, name = "~{basename}.lima.counts.txt" }
    call FF.FinalizeToFile as FinalizeDemultiplexReport { input: outdir = mdir, file = Demultiplex.report, name = "~{basename}.lima.report.txt" }
    call FF.FinalizeToFile as FinalizeDemultiplexSummary { input: outdir = mdir, file = Demultiplex.summary, name = "~{basename}.lima.summary.txt" }

    call FF.FinalizeToDir as FinalizeLimaSummary { input: outdir = fdir + "/summary/png", files = SummarizedDemuxReportPNG.report_files }
    call FF.FinalizeToDir as FinalizeLimaDetailedPNG { input: outdir = fdir + "/detailed/png", files = DetailedDemuxReportPNG.report_files }

    output {
        Array[File] demux_bams = FinalizeBam.gcs_path
        Array[File] demux_pbis = FinalizePbi.gcs_path

        File demux_counts = FinalizeDemultiplexCounts.gcs_path
        File demux_reports = FinalizeDemultiplexReport.gcs_path
        File demux_summary = FinalizeDemultiplexSummary.gcs_path

        String lima_summary_pngs = FinalizeLimaSummary.gcs_dir
        String lima_detailed_pngs = FinalizeLimaDetailedPNG.gcs_dir
    }
}
