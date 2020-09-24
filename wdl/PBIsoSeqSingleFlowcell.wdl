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
import "tasks/Finalize.wdl" as FF

workflow PBIsoSeqSingleFlowcell {
    input {
        String raw_reads_gcs_bucket
        String? sample_name

        File ref_fasta
        File ref_fasta_fai
        File ref_dict
        File ref_flat

        File barcode_file

        Int num_shards = 300

        String gcs_out_root_dir
    }

    parameter_meta {
        raw_reads_gcs_bucket: "GCS bucket holding subreads BAMs (and other related files) holding the sequences to be CCS-ed"
        sample_name:          "[optional] name of sample this FC is sequencing"

        ref_fasta:            "Reference fasta file"
        ref_fasta_fai:        "Index (.fai) for the reference fasta file"
        ref_dict:             "Sequence dictionary (.dict) for the reference fasta file"
        ref_flat:             "Gene predictions in refFlat format (https://genome.ucsc.edu/goldenpath/gbdDescriptions.html)"

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
        call Utils.ShardLongReads { input: unaligned_bam = subread_bam, unaligned_pbi = subread_pbi, num_shards = num_shards }

        # then perform correction on each of the shard
        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads }

            call AR.Minimap2 as AlignChunk {
                input:
                    reads      = [ CCS.consensus ],
                    ref_fasta  = ref_fasta,
                    RG         = RG,
                    map_preset = "splice"
            }
        }

        # merge the corrected per-shard BAM/report into one, corresponding to one raw input BAM
        call Utils.MergeBams as MergeChunks { input: bams = CCS.consensus, prefix = "~{SM}.~{ID}" }
        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report }

        call Utils.MergeBams as MergeAlignedChunks { input: bams = AlignChunk.aligned_bam, prefix = "~{SM}.~{ID}" }
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
            bam = ccs_bam,
            prefix = "~{SM[0]}.~{ID[0]}",
            barcode_file = barcode_file,
            isoseq = true,
            peek_guess = false
    }

    scatter (demux_bam in Demultiplex.demux_bams) {
        call PB.RefineTranscriptReads {
            input:
                bam = demux_bam,
                barcode_file = barcode_file,
                prefix = "~{SM[0]}.~{ID[0]}.flnc"
        }

        call PB.ClusterTranscripts {
            input:
                bam = RefineTranscriptReads.refined_bam,
                prefix = "~{SM[0]}.~{ID[0]}.clustered"
        }
    }

    scatter (p in zip(["refined", "clustered", "hq", "lq"],
                      [RefineTranscriptReads.refined_bam, ClusterTranscripts.clustered_bam, ClusterTranscripts.hq_bam, ClusterTranscripts.lq_bam])) {
        String RGA = "@RG\\tID:~{ID[0]}.~{p.left}\\tSM:~{SM[0]}\\tPL:~{PL[0]}\\tPU:~{PU[0]}\\tDT:~{DT[0]}"

        call AR.Minimap2 as AlignBAM {
            input:
                reads      = p.right,
                ref_fasta  = ref_fasta,
                RG         = RGA,
                map_preset = "splice",
                prefix     = "~{SM[0]}.~{ID[0]}.~{p.left}",
                runtime_attr_override = { "cpu_cores": 16 }
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
