version 1.0

##########################################################################################
## A workflow that performs demultiplexing on PacBio CLR reads from a single flow cell.
## The workflow demultiplexes CLR reads directly. A number of metrics and figures are
## produced along the way.
##########################################################################################

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/PBUtils.wdl" as PB
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/ShardUtils.wdl" as SU
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/Utils.wdl" as Utils
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/AlignReads.wdl" as AR
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/Finalize.wdl" as FF

workflow PBCLRDemultiplexOnlySingleFlowcell {
    input {
        String raw_reads_gcs_bucket
        String? sample_name
        File barcode_file
        String gcs_out_root_dir
    }

    parameter_meta {
        raw_reads_gcs_bucket: "GCS bucket holding subreads BAMs (and other related files) holding the sequences to be processed"
        sample_name:          "[optional] name of sample this FC is sequencing"
        barcode_file:         "GCS path to the fasta file that specifies the expected set of multiplexing barcodes"
        gcs_out_root_dir :    "GCS bucket to store the corrected/uncorrected reads and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = raw_reads_gcs_bucket }

    # double scatter: one FC may generate multiple raw BAMs, we perform another layer scatter on each of these BAMs
    scatter (subread_bam in FindBams.subread_bams) {
        call PB.GetRunInfo { input: subread_bam = subread_bam }

        String SM  = select_first([sample_name, GetRunInfo.run_info["SM"]])
        String PU  = GetRunInfo.run_info["PU"]
        String ID  = PU
        String DIR = SM + "." + ID

        # demultiplex CLR BAM
        call PB.Demultiplex { input: bam = subread_bam, ccs = false, guess = 0, prefix = "~{SM}.~{ID}", barcode_file = barcode_file }

        # make reports on demultiplexing
        call PB.MakeSummarizedDemultiplexingReport as SummarizedDemuxReportPNG { input: report = Demultiplex.report }
        call PB.MakeDetailedDemultiplexingReport as DetailedDemuxReportPNG { input: report = Demultiplex.report, type="png" }
        call PB.MakeDetailedDemultiplexingReport as DetailedDemuxReportPDF { input: report = Demultiplex.report, type="pdf" }
        call PB.MakePerBarcodeDemultiplexingReports as PerBarcodeDetailedDemuxReportPNG { input: report = Demultiplex.report, type="png" }
        call PB.MakePerBarcodeDemultiplexingReports as PerBarcodeDetailedDemuxReportPDF { input: report = Demultiplex.report, type="pdf" }

        ##########
        # store the results into designated bucket
        ##########

        call FF.FinalizeToDir as FinalizeDemuxReads {
            input:
                files = Demultiplex.demux_bams,
                outdir = outdir + "/" + DIR + "/reads"
        }

        call FF.FinalizeToDir as FinalizeLimaMetrics {
            input:
                files = [ Demultiplex.counts, Demultiplex.report, Demultiplex.summary ],
                outdir = outdir + "/" + DIR + "/metrics/lima"
        }

        call FF.FinalizeToDir as FinalizeLimaSummary {
            input:
                files = SummarizedDemuxReportPNG.report_files,
                outdir = outdir + "/" + DIR + "/figures/lima/summary/png"
        }

        call FF.FinalizeToDir as FinalizeLimaDetailedPNG {
            input:
                files = DetailedDemuxReportPNG.report_files,
                outdir = outdir + "/" + DIR + "/figures/lima/detailed/png"
        }

        call FF.FinalizeToDir as FinalizeLimaDetailedPDF {
            input:
                files = DetailedDemuxReportPDF.report_files,
                outdir = outdir + "/" + DIR + "/figures/lima/detailed/pdf"
        }
    }
}
