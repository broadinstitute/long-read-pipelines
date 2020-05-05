version 1.0

##########################################################################################
## A workflow that performs CCS correction and IsoSeq processing on PacBio HiFi reads from
## a single flow cell. The workflow shards the subreads into clusters and performs CCS in
## parallel on each cluster.  Error-corrected reads are then processed with PacBio's
## IsoSeq software.  A number of metrics and figures are produced along the way.
##########################################################################################

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/Tama.wdl" as TAMA
import "tasks/Finalize.wdl" as FF

workflow PBIsoSeqSingleFlowcell {
    input {
        String raw_reads_gcs_bucket
        String? sample_name

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File barcode_file = "gs://broad-dsde-methods-long-reads/resources/multiplexing/isoseq_primers_broad.fasta"

        Int num_shards = 300

        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/PBIsoSeqSingleFlowcell"
    }

    parameter_meta {
        raw_reads_gcs_bucket: "GCS bucket holding subreads BAMs (and other related files) holding the sequences to be CCS-ed"
        sample_name:          "[optional] name of sample this FC is sequencing"

        ref_fasta:            "Reference fasta file"
        ref_fasta_fai:        "Index (.fai) for the reference fasta file"
        ref_dict:             "Sequence dictionary (.dict) for the reference fasta file"

        barcode_file:         "GCS path to the fasta file that specifies the expected set of multiplexing barcodes"

        num_shards:           "[default-valued] number of sharded BAMs to create (tune for performance)"

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
        String RG = "@RG\\tID:~{ID}\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        # break one raw BAM into fixed number of shards
        File subread_pbi = sub(subread_bam, ".bam$", ".bam.pbi")
        call PB.ShardLongReads {
            input:
                unaligned_bam = subread_bam,
                unaligned_pbi = subread_pbi,
                num_shards = num_shards
        }

        # then perform correction on each of the shard
        RuntimeAttr ccs_runtime_attrs = object {
            mem_gb: 16,
            preemptible_tries: 0
        }
        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads, runtime_attr_override = ccs_runtime_attrs }
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
    call PB.Demultiplex {
        input:
            bam          = ccs_bam,
            prefix       = "~{SM[0]}.~{ID[0]}",
            barcode_file = barcode_file,
            isoseq       = true,
            peek_guess   = false
    }

    scatter (demux_bam in Demultiplex.demux_bams) {
        call PB.RefineTranscriptReads {
            input:
                bam          = demux_bam,
                barcode_file = barcode_file,
                prefix       = "~{SM[0]}.~{ID[0]}.flnc"
        }

        call PB.ClusterTranscripts {
            input:
                bam          = RefineTranscriptReads.refined_bam,
                prefix       = "~{SM[0]}.~{ID[0]}.clustered"
        }

        call PB.Align as AlignTranscripts {
            input:
                bam          = ClusterTranscripts.clustered_bam,
                ref_fasta    = ref_fasta,
                sample_name  = SM[0],
                map_preset   = "ISOSEQ",
                prefix       = "~{SM[0]}.~{ID[0]}",
                runtime_attr_override = { "cpu_cores": 32 }
        }

        call PB.CollapseTranscripts {
            input:
                bam          = AlignTranscripts.aligned_bam,
                prefix       = "~{SM[0]}.~{ID[0]}.collapsed"
        }

        call TAMA.CollapseIsoforms {
            input:
                bam          = AlignTranscripts.aligned_bam,
                ref_fasta    = ref_fasta,
                prefix       = "~{SM[0]}.~{ID[0]}.collapsed"
        }
    }

    ##########
    # store the results into designated bucket
    ##########

#    call FF.FinalizeToDir as FinalizeMergedRuns {
#        input:
#            files = [ ccs_bam ],
#            outdir = outdir + "/" + DIR[0] + "/alignments"
#    }
#
#    call FF.FinalizeToDir as FinalizeCCSMetrics {
#        input:
#            files = [ ccs_report ],
#            outdir = outdir + "/" + DIR[0] + "/metrics/ccs_metrics"
#    }
}
