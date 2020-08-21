version 1.0

##########################################################################################
## A workflow that performs CCS correction, demultiplexing, and variant calling on PacBio
## HiFi reads from a single flow cell. The workflow shards the subreads into clusters and
## performs CCS in parallel on each cluster.  Error-corrected reads are then demultiplexed,
## and each demultiplexed dataset is independently variant-called.  A number of metrics
## and figures are produced along the way.
##########################################################################################

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.29/wdl/tasks/PBUtils.wdl" as PB
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.29/wdl/tasks/ShardUtils.wdl" as SU
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.29/wdl/tasks/Utils.wdl" as Utils
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.29/wdl/tasks/AlignReads.wdl" as AR
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.29/wdl/tasks/AlignedMetrics.wdl" as AM
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.29/wdl/tasks/Finalize.wdl" as FF
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.29/wdl/tasks/CallSVs.wdl" as SV
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.29/wdl/tasks/CallSmallVariants.wdl" as SMV

workflow PBCCSDemultiplexWholeGenomeSingleFlowcell {
    input {
        String raw_reads_gcs_bucket
        String? sample_name

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File tandem_repeat_bed
        File ref_flat
        File dbsnp_vcf
        File dbsnp_tbi

        String mt_chr_name
        File metrics_locus

        Int? num_reads_per_split = 200000
        File barcode_file

        String gcs_out_root_dir
    }

    parameter_meta {
        raw_reads_gcs_bucket: "GCS bucket holding subreads BAMs (and other related files) holding the sequences to be CCS-ed"
        sample_name:          "[optional] name of sample this FC is sequencing"

        ref_fasta:            "Reference fasta file"
        ref_fasta_fai:        "Index (.fai) for the reference fasta file"
        ref_dict:             "Sequence dictionary (.dict) for the reference fasta file"

        tandem_repeat_bed:    "BED file specifying the location of tandem repeats in the reference"
        ref_flat:             "Gene predictions in refFlat format (https://genome.ucsc.edu/goldenpath/gbdDescriptions.html)"
        dbsnp_vcf:            "dbSNP vcf"
        dbsnp_tbi:            "Index (.tbi) for dbSNP vcf"

        mt_chr_name:          "Contig name for the mitochondrial sequence in the reference"
        metrics_locus:        "Loci over which some summary metrics should be computed"

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
        String PL  = "PACBIO"
        String PU  = GetRunInfo.run_info["PU"]
        String DT  = GetRunInfo.run_info["DT"]
        String ID  = PU
        String DS  = GetRunInfo.run_info["DS"]
        String DIR = SM + "." + ID
        String RG  = "@RG\\tID:~{ID}\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        # shard one raw BAM into fixed chunk size (num_reads_per_split)
        call Utils.ShardLongReads { input: unmapped_files = [ subread_bam ], num_reads_per_split = 200000 }

        # then perform correction on each of the shard
        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads }
        }

        # merge the corrected per-shard BAM/report into one, corresponding to one raw input BAM
        call Utils.MergeBams as MergeChunks { input: bams = CCS.consensus }
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
    call PB.MakePerBarcodeDemultiplexingReports as PerBarcodeDetailedDemuxReportPNG { input: report = Demultiplex.report, type="png" }
    call PB.MakePerBarcodeDemultiplexingReports as PerBarcodeDetailedDemuxReportPDF { input: report = Demultiplex.report, type="pdf" }

    # scatter on each demultiplexed BAM file
    scatter (demux_bam in Demultiplex.demux_bams) {
        String BC   = sub(basename(demux_bam, ".bam"), "~{SM[0]}.~{ID[0]}.", "")
        String BCID = ID[0] + "." + BC
        String BCRG = "@RG\\tID:~{BCID}\\tSM:~{SM[0]}\\tPL:~{PL[0]}\\tPU:~{PU[0]}\\tDT:~{DT[0]}"

        # align reads to reference
        call AR.Minimap2 as AlignBarcode {
            input:
                reads      = [ demux_bam ],
                ref_fasta  = ref_fasta,
                RG         = BCRG,
                map_preset = "asm20",
                prefix     = BCID + ".aligned"
        }

        # compute alignment metrics
        call AM.AlignedMetrics as PerBarcodeMetrics {
            input:
                aligned_bam    = AlignBarcode.aligned_bam,
                aligned_bai    = AlignBarcode.aligned_bai,
                ref_fasta      = ref_fasta,
                ref_dict       = ref_dict,
                ref_flat       = ref_flat,
                dbsnp_vcf      = dbsnp_vcf,
                dbsnp_tbi      = dbsnp_tbi,
                metrics_locus  = metrics_locus,
                per            = "flowcell",
                type           = "barcode",
                label          = BC,
                gcs_output_dir = outdir + "/" + DIR[0]
        }

        # create a BED file that indicates where the BAM file has coverage
        call Utils.BamToBed { input: bam = AlignBarcode.aligned_bam, prefix = BCID }

        # call SVs
        call SV.CallSVs as CallSVs {
            input:
                bam               = AlignBarcode.aligned_bam,
                bai               = AlignBarcode.aligned_bai,

                ref_fasta         = ref_fasta,
                ref_fasta_fai     = ref_fasta_fai,
                tandem_repeat_bed = tandem_repeat_bed
        }

        # call SNVs and small indels
        call SMV.CallSmallVariants as CallSmallVariants {
            input:
                bam               = AlignBarcode.aligned_bam,
                bai               = AlignBarcode.aligned_bai,

                ref_fasta         = ref_fasta,
                ref_fasta_fai     = ref_fasta_fai,
                ref_dict          = ref_dict
        }

        ##########
        # store the demultiplexing results into designated bucket
        ##########

        call FF.FinalizeToDir as FinalizeDemuxAlignedReads {
            input:
                files  = [ AlignBarcode.aligned_bam, AlignBarcode.aligned_bai, BamToBed.bed ],
                outdir = outdir + "/" + DIR[0] + "/alignments/" + BC
        }

        call FF.FinalizeToDir as FinalizeSVs {
            input:
                files = [ CallSVs.pbsv_vcf, CallSVs.sniffles_vcf, CallSVs.svim_vcf ],
                outdir = outdir + "/" + DIR[0] + "/variants/" + BC
        }

        call FF.FinalizeToDir as FinalizeSmallVariants {
            input:
                files = [ CallSmallVariants.longshot_vcf, CallSmallVariants.longshot_tbi ],
                outdir = outdir + "/" + DIR[0] + "/variants/" + BC
        }
    }

    # merge demultiplexed BAMs into a single BAM (one readgroup per file)
    call Utils.MergeBams as MergeBCRuns { input: bams = AlignBarcode.aligned_bam, prefix = "barcodes" }

    ##########
    # store the results into designated bucket
    ##########

    call FF.FinalizeToDir as FinalizeDemuxCombinedReads {
        input:
            files = [ MergeBCRuns.merged_bam, MergeBCRuns.merged_bai ],
            outdir = outdir + "/" + DIR[0] + "/reads"
    }

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

    call FF.FinalizeToDir as FinalizeLimaPerBarcodeDetailedPNG {
        input:
            files = PerBarcodeDetailedDemuxReportPNG.report_files,
            outdir = outdir + "/" + DIR[0] + "/figures/lima/per_barcode/png"
    }

    call FF.FinalizeToDir as FinalizeLimaPerBarcodeDetailedPDF {
        input:
            files = PerBarcodeDetailedDemuxReportPDF.report_files,
            outdir = outdir + "/" + DIR[0] + "/figures/lima/per_barcode/pdf"
    }
}
