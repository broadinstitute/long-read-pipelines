version 1.0

import "../../../tasks/Utility/PBUtils.wdl" as PB
import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow PBCCSIsoSeq {

    meta {
        description: "A workflow that performs CCS correction and IsoSeq processing on PacBio HiFi reads from a single flow cell. The workflow shards the subreads into clusters and performs CCS in parallel on each cluster. Error-corrected reads are then processed with PacBio's IsoSeq software. A number of metrics and figures are produced along the way."
    }
    parameter_meta {
        ccs_bams:         "GCS path to CCS BAM files"
        ccs_pbis:         "GCS path to CCS BAM .pbi indices"

        ref_map_file:     "table indicating reference sequence and auxillary file locations"
        participant_name: "name of the participant from whom these samples were obtained"
        barcode_file:     "GCS path to the fasta file that specifies the expected set of multiplexing barcodes"

        gcs_out_root_dir: "GCS bucket to store the corrected/uncorrected reads, variants, and metrics files"
    }

    input {
        Array[File] ccs_bams
        Array[File] ccs_pbis

        File ref_map_file
        String participant_name
        File barcode_file

        Boolean drop_per_base_N_pulse_tags = true

        String gcs_out_root_dir
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCCSIsoSeq/~{participant_name}"

    # gather across (potential multiple) input CCS BAMs
    if (length(ccs_bams) > 1) {
        call Utils.MergeBams as MergeAllReads { input: bams = ccs_bams, prefix = participant_name }
        call PB.PBIndex as IndexCCSUnalignedReads { input: bam = MergeAllReads.merged_bam }
    }

    File bam = select_first([MergeAllReads.merged_bam, ccs_bams[0]])
    File pbi = select_first([IndexCCSUnalignedReads.pbi, ccs_pbis[0]])

    # demultiplex CCS-ed BAM
    call PB.Demultiplex {
        input:
            bam = bam,
            prefix = participant_name,
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
        String BC = sub(basename(demux_bam, ".bam"), "~{participant_name}.corrected", "")

        call PB.RefineTranscriptReads {
            input:
                bam          = demux_bam,
                barcode_file = barcode_file,
                prefix       = "~{participant_name}.~{BC}.flnc"
        }

        call PB.ClusterTranscripts {
            input:
                bam          = RefineTranscriptReads.refined_bam,
                prefix       = "~{participant_name}.~{BC}.clustered"
        }

        call PB.Align as AlignTranscripts {
            input:
                bam          = ClusterTranscripts.clustered_bam,
                ref_fasta    = ref_map['fasta'],
                sample_name  = participant_name,
                drop_per_base_N_pulse_tags = drop_per_base_N_pulse_tags,
                map_preset   = "ISOSEQ",
                prefix       = "~{participant_name}.~{BC}",
                runtime_attr_override = { "cpu_cores": 32 }
        }

        # create a BED file that indicates where the BAM file has coverage
        call Utils.BamToBed { input: bam = AlignTranscripts.aligned_bam, prefix = BC }

        call PB.CollapseTranscripts {
            input:
                bam          = AlignTranscripts.aligned_bam,
                prefix       = "~{participant_name}.~{BC}.collapsed"
        }

        ##########
        # store the demultiplexing results into designated bucket
        ##########

        String adir = outdir + "/" + BC + "/alignments"
        String tdir = outdir + "/" + BC + "/transcripts"

        call FF.FinalizeToFile as FinalizeAlignedTranscriptsBam { input: outdir = adir, file = AlignTranscripts.aligned_bam }
        call FF.FinalizeToFile as FinalizeAlignedTranscriptsBai { input: outdir = adir, file = AlignTranscripts.aligned_bai }
        call FF.FinalizeToFile as FinalizeAlignedTranscriptsBed { input: outdir = adir, file = BamToBed.bed }
        call FF.FinalizeToFile as FinalizeCollapsedTranscripts { input: outdir = tdir, file = CollapseTranscripts.gff }
    }

    # merge demultiplexed BAMs into a single BAM (one readgroup per file)
    call Utils.MergeBams as MergeBarcodeBams { input: bams = AlignTranscripts.aligned_bam, prefix = "barcodes" }

    call PB.PBIndex as IndexAlignedReads { input: bam = MergeBarcodeBams.merged_bam }

    # Finalize
    String rdir = outdir + "/reads"
    String bdir = outdir + "/alignments/all_barcodes"
    String mdir = outdir + "/metrics/combined/lima"
    String fdir = outdir + "/figures"

    call FF.FinalizeToFile as FinalizeBam { input: outdir = rdir, file = bam, name = "~{participant_name}.bam" }
    call FF.FinalizeToFile as FinalizePbi { input: outdir = rdir, file = pbi, name = "~{participant_name}.bam.pbi" }

    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = bdir, file = MergeBarcodeBams.merged_bam, name = "~{participant_name}.all_barcodes.bam" }
    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = bdir, file = MergeBarcodeBams.merged_bai, name = "~{participant_name}.all_barcodes.bam.bai"  }
    call FF.FinalizeToFile as FinalizeAlignedPbi { input: outdir = bdir, file = IndexAlignedReads.pbi, name = "~{participant_name}.all_barcodes.bam.pbi"  }

    call FF.FinalizeToFile as FinalizeDemultiplexCounts { input: outdir = mdir, file = Demultiplex.counts, name = "~{participant_name}.lima.counts.txt" }
    call FF.FinalizeToFile as FinalizeDemultiplexReport { input: outdir = mdir, file = Demultiplex.report, name = "~{participant_name}.lima.report.txt" }
    call FF.FinalizeToFile as FinalizeDemultiplexSummary { input: outdir = mdir, file = Demultiplex.summary, name = "~{participant_name}.lima.summary.txt" }

    call FF.FinalizeToDir as FinalizeLimaSummary { input: outdir = fdir + "/summary/png", files = SummarizedDemuxReportPNG.report_files }
    call FF.FinalizeToDir as FinalizeLimaDetailedPNG { input: outdir = fdir + "/detailed/png", files = DetailedDemuxReportPNG.report_files }

    output {
        File ccs_bam = FinalizeBam.gcs_path
        File ccs_pbi = FinalizePbi.gcs_path

        File aligned_bam = FinalizeAlignedBam.gcs_path
        File aligned_bai = FinalizeAlignedBai.gcs_path
        File aligned_pbi = FinalizeAlignedPbi.gcs_path

        File demux_counts = FinalizeDemultiplexCounts.gcs_path
        File demux_reports = FinalizeDemultiplexReport.gcs_path
        File demux_summary = FinalizeDemultiplexSummary.gcs_path
    }
}
