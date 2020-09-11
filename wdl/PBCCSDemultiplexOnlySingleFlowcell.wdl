version 1.0

##########################################################################################
## A workflow that performs CCS correction and demultiplexing on PacBio HiFi reads from a
## single flow cell. The workflow shards the subreads into clusters and performs CCS in
## parallel on each cluster.  Error-corrected reads are then demultiplexed.  A number of
## metrics and figures are produced along the way.
##########################################################################################

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.30/wdl/tasks/PBUtils.wdl" as PB
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.30/wdl/tasks/ShardUtils.wdl" as SU
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.30/wdl/tasks/Utils.wdl" as Utils
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.30/wdl/tasks/Finalize.wdl" as FF

workflow PBCCSDemultiplexOnlySingleFlowcell {
    input {
        String raw_reads_gcs_bucket

        String? sample_name
        Int num_reads_per_split = 2000000

        File barcode_file

        String gcs_out_root_dir
    }

    parameter_meta {
        raw_reads_gcs_bucket: "GCS bucket holding subreads BAMs (and other related files) holding the sequences to be CCS-ed"
        sample_name:          "[optional] name of sample this FC is sequencing"
        num_reads_per_split:  "[default-valued] number of subreads each sharded BAM contains (tune for performance)"
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

        # shard one raw BAM into fixed chunk size (num_reads_per_split)
        call Utils.ShardLongReads { input: unmapped_files = [ subread_bam ], num_reads_per_split = num_reads_per_split }

        # then perform correction on each of the shard
        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads }
        }

        # merge the corrected per-shard BAM/report into one, corresponding to one raw input BAM
        call Utils.MergeBams as MergeChunks { input: bams = CCS.consensus, prefix = "~{SM}.~{ID}" }
        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report }
    }

    # gather across (potential multiple) input raw BAMs
    if (length(FindBams.subread_bams) > 1) {
        call Utils.MergeBams as MergeRuns { input: bams = MergeChunks.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }
        call PB.MergeCCSReports as MergeAllCCSReports { input: reports = MergeCCSReports.report }
    }

    File ccs_bam = select_first([ MergeRuns.merged_bam, MergeChunks.merged_bam[0] ])
    File ccs_report = select_first([ MergeAllCCSReports.report, MergeCCSReports.report[0] ])

    # demultiplex CCS-ed BAM
    call PB.Demultiplex { input: bam = ccs_bam, prefix = "~{SM[0]}.~{ID[0]}", barcode_file = barcode_file }

    # make reports on demultiplexing
    call PB.MakeSummarizedDemultiplexingReport as SummarizedDemuxReportPNG { input: report = Demultiplex.report }
    call PB.MakeDetailedDemultiplexingReport as DetailedDemuxReportPNG { input: report = Demultiplex.report, type="png" }
    call PB.MakeDetailedDemultiplexingReport as DetailedDemuxReportPDF { input: report = Demultiplex.report, type="pdf" }

    ##########
    # store the results into designated bucket
    ##########

    call FF.FinalizeToDir as FinalizeDemuxReads {
        input:
            files = Demultiplex.demux_bams,
            outdir = outdir + "/" + DIR[0] + "/reads"
    }

    call FF.FinalizeToDir as FinalizeCCSMetrics {
        input:
            files = [ ccs_report ],
            outdir = outdir + "/" + DIR[0] + "/metrics/ccs"
    }

    call FF.FinalizeToDir as FinalizeLimaMetrics {
        input:
            files = [ Demultiplex.counts, Demultiplex.report, Demultiplex.summary ],
            outdir = outdir + "/" + DIR[0] + "/metrics/lima"
    }

    call FF.FinalizeToDir as FinalizeLimaSummary {
        input:
            files = SummarizedDemuxReportPNG.report_files,
            outdir = outdir + "/" + DIR[0] + "/figures/lima/summary/png"
    }

    call FF.FinalizeToDir as FinalizeLimaDetailedPNG {
        input:
            files = DetailedDemuxReportPNG.report_files,
            outdir = outdir + "/" + DIR[0] + "/figures/lima/detailed/png"
    }

    call FF.FinalizeToDir as FinalizeLimaDetailedPDF {
        input:
            files = DetailedDemuxReportPDF.report_files,
            outdir = outdir + "/" + DIR[0] + "/figures/lima/detailed/pdf"
    }
}
