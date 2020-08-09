version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/ShardUtils.wdl" as SU
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/Finalize.wdl" as FF
import "tasks/CallSVs.wdl" as SV
import "tasks/CallSmallVariants.wdl" as SMV

workflow PBCCSDemultiplexWholeGenomeSingleFlowcell {
    input {
        String gcs_input_dir
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

        File barcode_file

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = gcs_input_dir }

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

        call Utils.ShardLongReads { input: unmapped_files = [ subread_bam ], num_reads_per_split = 200000 }

        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads }
        }

        call Utils.MergeBams as MergeChunks { input: bams = CCS.consensus }
        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report }
    }

    if (length(FindBams.subread_bams) > 1) {
        call Utils.MergeBams as MergeRuns { input: bams = MergeChunks.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }
        call PB.MergeCCSReports as MergeAllCCSReports { input: reports = MergeCCSReports.report }
    }

    File ccs_bam = select_first([ MergeRuns.merged_bam, MergeChunks.merged_bam[0] ])
    File ccs_report = select_first([ MergeAllCCSReports.report, MergeCCSReports.report[0] ])

    call PB.Demultiplex { input: bam = ccs_bam, prefix = "~{SM[0]}.~{ID[0]}", barcode_file = barcode_file }

    call PB.MakeSummarizedDemultiplexingReport as SummarizedDemuxReportPNG { input: report = Demultiplex.report }
    call PB.MakeDetailedDemultiplexingReport as DetailedDemuxReportPNG { input: report = Demultiplex.report, type="png" }
    call PB.MakeDetailedDemultiplexingReport as DetailedDemuxReportPDF { input: report = Demultiplex.report, type="pdf" }
    call PB.MakePerBarcodeDemultiplexingReports as PerBarcodeDetailedDemuxReportPNG { input: report = Demultiplex.report, type="png" }
    call PB.MakePerBarcodeDemultiplexingReports as PerBarcodeDetailedDemuxReportPDF { input: report = Demultiplex.report, type="pdf" }

    scatter (demux_bam in Demultiplex.demux_bams) {
        String BC   = sub(basename(demux_bam, ".bam"), "~{SM[0]}.~{ID[0]}.", "")
        String BCID = ID[0] + "." + BC
        String BCRG = "@RG\\tID:~{BCID}\\tSM:~{SM[0]}\\tPL:~{PL[0]}\\tPU:~{PU[0]}\\tDT:~{DT[0]}"

        call AR.Minimap2 as AlignBarcode {
            input:
                reads      = [ demux_bam ],
                ref_fasta  = ref_fasta,
                RG         = BCRG,
                map_preset = "asm20",
                prefix     = BCID + ".aligned"
        }

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

        call Utils.BamToBed { input: bam = AlignBarcode.aligned_bam, prefix = BCID }

        call SV.CallSVs as CallSVs {
            input:
                bam               = AlignBarcode.aligned_bam,
                bai               = AlignBarcode.aligned_bai,

                ref_fasta         = ref_fasta,
                ref_fasta_fai     = ref_fasta_fai,
                tandem_repeat_bed = tandem_repeat_bed
        }

        call SMV.CallSmallVariants as CallSmallVariants {
            input:
                bam               = AlignBarcode.aligned_bam,
                bai               = AlignBarcode.aligned_bai,

                ref_fasta         = ref_fasta,
                ref_fasta_fai     = ref_fasta_fai,
                ref_dict          = ref_dict
        }

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

    call Utils.MergeBams as MergeBCRuns { input: bams = AlignBarcode.aligned_bam, prefix = "barcodes" }

    ##########
    # Finalize
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
            files = [ Demultiplex.counts, Demultiplex.guess, Demultiplex.report, Demultiplex.summary ],
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
