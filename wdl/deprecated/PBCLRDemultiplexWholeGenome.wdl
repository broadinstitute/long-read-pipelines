version 1.0

##########################################################################################
## A workflow that performs demultiplexing on PacBio CLR reads from a single flow cell.
## The workflow demultiplexes CLR reads directly. A number of metrics and figures are
## produced along the way.
##########################################################################################

import "tasks/Utility/PBUtils.wdl" as DeprecatedPB
import "../tasks/Utility/PBUtils.wdl" as PB
import "../tasks/Utility/Utils.wdl" as Utils
import "../tasks/Alignment/AlignReads.wdl" as AR
import "../tasks/QC/AlignedMetrics.wdl" as AM
import "../tasks/Utility/Finalize.wdl" as FF

import "tasks/CallVariantsPBCLR.wdl" as VAR

workflow PBCLRDemultiplexWholeGenome {
    input {
        Array[File] bams
        File ref_map_file

        String participant_name
        File barcode_file
        Int num_shards = 300

        Boolean drop_per_base_N_pulse_tags = true

        String gcs_out_root_dir
    }

    parameter_meta {
        bams:             "GCS path to raw subreads or CCS data"
        ref_map_file:     "table indicating reference sequence and auxillary file locations"

        participant_name: "name of the participant from whom these samples were obtained"
        barcode_file:     "GCS path to the fasta file that specifies the expected set of multiplexing barcodes"
        num_shards:       "[default-valued] number of sharded BAMs to create (tune for performance)"

        gcs_out_root_dir: "GCS bucket to store the corrected/uncorrected reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCLRDemultiplexWholeGenome/" + participant_name

    # gather across (potential multiple) input raw BAMs
    if (length(bams) > 1) {
        call Utils.MergeBams as MergeAllUncorrected { input: bams = bams, prefix = "~{participant_name}.uncorrected" }
    }

    File uncorrected_bam = select_first([ MergeAllUncorrected.merged_bam, bams[0] ])

    # demultiplex CLR BAM
    call PB.Demultiplex {
        input:
            bam = uncorrected_bam,
            prefix = participant_name,
            barcode_file = barcode_file,
            ccs = false,
            guess = 75,
            guess_min_count = 1,
            dump_removed = true,
            split_bam_named = true
    }

    # make reports on demultiplexing
    call PB.MakeSummarizedDemultiplexingReport as SummarizedDemuxReportPNG { input: report = Demultiplex.report }
    call PB.MakeDetailedDemultiplexingReport as DetailedDemuxReportPNG { input: report = Demultiplex.report, type="png" }
    call PB.MakeDetailedDemultiplexingReport as DetailedDemuxReportPDF { input: report = Demultiplex.report, type="pdf" }
    call DeprecatedPB.MakePerBarcodeDemultiplexingReports as PerBarcodeDetailedDemuxReportPNG { input: report = Demultiplex.report, type="png" }
    call DeprecatedPB.MakePerBarcodeDemultiplexingReports as PerBarcodeDetailedDemuxReportPDF { input: report = Demultiplex.report, type="pdf" }

    # scatter on each demultiplexed BAM file
    scatter (demux_bam in Demultiplex.demux_bams) {
        String BC = sub(basename(demux_bam, ".bam"), "~{participant_name}.uncorrected", "")

        # align reads to reference
        call PB.Align as AlignBarcode {
            input:
                bam         = demux_bam,
                ref_fasta   = ref_map['fasta'],
                sample_name = participant_name,
                map_preset  = "CLR",
                drop_per_base_N_pulse_tags = drop_per_base_N_pulse_tags
        }

        # compute alignment metrics
        call AM.AlignedMetrics as PerBarcodeMetrics {
            input:
                aligned_bam    = AlignBarcode.aligned_bam,
                aligned_bai    = AlignBarcode.aligned_bai,
                ref_fasta      = ref_map['fasta'],
                ref_dict       = ref_map['dict'],
                gcs_output_dir = outdir + "/" + BC + "/metrics/per_barcode"
        }

        # create a BED file that indicates where the BAM file has coverage
        call Utils.BamToBed { input: bam = AlignBarcode.aligned_bam, prefix = BC }

        call VAR.CallVariants {
            input:
                bam               = AlignBarcode.aligned_bam,
                bai               = AlignBarcode.aligned_bai,

                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                tandem_repeat_bed = ref_map['tandem_repeat_bed'],
        }

        ##########
        # store the demultiplexing results into designated bucket
        ##########

        call FF.FinalizeToDir as FinalizeDemuxAlignedReads {
            input:
                files  = [ AlignBarcode.aligned_bam, AlignBarcode.aligned_bai, BamToBed.bed ],
                outdir = outdir + "/" + BC + "/alignments"
        }

    }

    # merge demultiplexed BAMs into a single BAM (one readgroup per file)
    call Utils.MergeBams as MergeBarcodeBams { input: bams = AlignBarcode.aligned_bam, prefix = "barcodes" }

    ##########
    # store the results into designated bucket
    ##########

    call FF.FinalizeToDir as FinalizeDemuxCombinedReads {
        input:
            files = [ MergeBarcodeBams.merged_bam, MergeBarcodeBams.merged_bai ],
            outdir = outdir + "/combined/alignments"
    }

    call FF.FinalizeToDir as FinalizeLimaMetrics {
        input:
            files = [ Demultiplex.counts, Demultiplex.report, Demultiplex.summary ],
            outdir = outdir + "/metrics/combined/" + participant_name + "/lima"
    }

    call FF.FinalizeToDir as FinalizeLimaSummary {
        input:
            files = SummarizedDemuxReportPNG.report_files,
            outdir = outdir + "/figures/combined/" + participant_name + "/lima/summary/png"
    }

    call FF.FinalizeToDir as FinalizeLimaDetailedPNG {
        input:
            files = DetailedDemuxReportPNG.report_files,
            outdir = outdir + "/figures/combined/" + participant_name + "/lima/detailed/png"
    }

    call FF.FinalizeToDir as FinalizeLimaDetailedPDF {
        input:
            files = DetailedDemuxReportPDF.report_files,
            outdir = outdir + "/figures/combined/" + participant_name + "/lima/detailed/pdf"
    }

    call FF.FinalizeToDir as FinalizeLimaPerBarcodeDetailedPNG {
        input:
            files = PerBarcodeDetailedDemuxReportPNG.report_files,
            outdir = outdir + "/figures/per_barcode/lima/png"
    }

    call FF.FinalizeToDir as FinalizeLimaPerBarcodeDetailedPDF {
        input:
            files = PerBarcodeDetailedDemuxReportPDF.report_files,
            outdir = outdir + "/figures/per_barcode/lima/pdf"
    }
}
