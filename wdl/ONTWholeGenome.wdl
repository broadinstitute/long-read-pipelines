version 1.0

import "tasks/ONTUtils.wdl" as ONT
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/CallSVs.wdl" as SV
import "tasks/Figures.wdl" as FIG
import "tasks/Finalize.wdl" as FF
import "tasks/CallSmallVariantsONT.wdl" as SMV
import "tasks/Methylation.wdl" as Meth

workflow ONTWholeGenome {
    input {
        Array[File] final_summaries
        Array[File] sequencing_summaries
        File ref_map_file

        String participant_name
        Int num_shards = 50

        String gcs_out_root_dir
    }

    parameter_meta {
        final_summaries:      "GCS path to '*final_summary*.txt*' files for basecalled fastq files"
        sequencing_summaries: "GCS path to '*sequencing_summary*.txt*' files for basecalled fastq files"
        ref_map_file:         "table indicating reference sequence and auxillary file locations"

        participant_name:     "name of the participant from whom these samples were obtained"
        num_shards:           "[default-valued] number of shards into which fastq files should be batched"

        gcs_out_root_dir:     "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTWholeGenome/" + participant_name

    scatter (p in zip(final_summaries, sequencing_summaries)) {
        File final_summary = p.left
        File sequencing_summary = p.right

        call ONT.GetRunInfo { input: summary_file = final_summary }
        call ONT.ListFiles as ListFast5s { input: summary_file = final_summary, suffix = "fast5" }
        call ONT.ListFiles as ListFastqs { input: summary_file = final_summary, suffix = "fastq" }

        String SM  = participant_name
        String PL  = "ONT"
        String PU  = GetRunInfo.run_info["instrument"]
        String DT  = GetRunInfo.run_info["started"]
        String ID  = GetRunInfo.run_info["flow_cell_id"] + "." + GetRunInfo.run_info["position"]
        String DIR = GetRunInfo.run_info["protocol_group_id"] + "." + SM + "." + ID
        String SID = ID + "." + sub(GetRunInfo.run_info["protocol_run_id"], "-.*", "")
        String RG = "@RG\\tID:~{SID}\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        call ONT.PartitionManifest as PartitionFastqManifest { input: manifest = ListFastqs.manifest, N = num_shards }

        scatter (manifest_chunk in PartitionFastqManifest.manifest_chunks) {
            call AR.Minimap2 as AlignReads {
                input:
                    reads      = read_lines(manifest_chunk),
                    ref_fasta  = ref_map['fasta'],
                    RG         = RG,
                    map_preset = "map-ont"
            }
        }

        call Utils.MergeBams as MergeReads { input: bams = AlignReads.aligned_bam }

#        call Meth.Methylation {
#            input:
#                fast5s = read_lines(ListFast5s.manifest),
#                fastqs = read_lines(ListFastqs.manifest),
#                sequencing_summary = sequencing_summary,
#
#                bam = MergeReads.merged_bam,
#                bai = MergeReads.merged_bai,
#
#                ref_fasta = ref_map['fasta'],
#                ref_fai   = ref_map['fai'],
#
#                prefix = "~{SM}.~{ID}.methylation"
#        }

        call AM.AlignedMetrics as PerFlowcellMetrics {
            input:
                aligned_bam    = MergeReads.merged_bam,
                aligned_bai    = MergeReads.merged_bai,
                ref_fasta      = ref_map['fasta'],
                ref_dict       = ref_map['dict'],
                gcs_output_dir = outdir + "/metrics/per_flowcell/" + SID
        }

        call FIG.Figures as PerFlowcellFigures {
            input:
                summary_files  = [ sequencing_summary ],

                gcs_output_dir = outdir + "/metrics/per_flowcell/" + SID
        }
    }

    if (length(final_summaries) > 1) {
        call Utils.MergeBams as MergeAllReads { input: bams = MergeReads.merged_bam, prefix = "~{participant_name}" }
    }

    File bam = select_first([ MergeAllReads.merged_bam, MergeReads.merged_bam[0] ])
    File bai = select_first([ MergeAllReads.merged_bai, MergeReads.merged_bai[0] ])

    call AM.AlignedMetrics as PerFlowcellRunMetrics {
        input:
            aligned_bam    = bam,
            aligned_bai    = bai,
            ref_fasta      = ref_map['fasta'],
            ref_dict       = ref_map['dict'],
            gcs_output_dir = outdir + "/metrics/combined/" + participant_name
    }

    call FIG.Figures as PerSampleFigures {
        input:
            summary_files  = sequencing_summaries,

            gcs_output_dir = outdir + "/metrics/combined/" + participant_name
    }

#    call SV.CallSVs as CallSVs {
#        input:
#            bam               = bam,
#            bai               = bai,
#
#            ref_fasta         = ref_map['fasta'],
#            ref_fasta_fai     = ref_map['fai'],
#            tandem_repeat_bed = ref_map['tandem_repeat_bed'],
#
#            preset            = "ont"
#    }

    call SMV.CallSmallVariants as CallSmallVariants {
        input:
            bam               = bam,
            bai               = bai,

            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
            ref_dict          = ref_map['dict'],
    }

#    if (length(Methylation.freq_tsv) > 1) {
#        call Meth.FreqMerge { input: freq_tsvs = Methylation.freq_tsv, prefix = "~{SM[0]}.~{ID[0]}.methylation" }
#    }
#    File freq_tsv = select_first([ FreqMerge.freq_tsv , Methylation.freq_tsv[0] ])

    ##########
    # Finalize
    ##########

#    call FF.FinalizeToDir as FinalizeSVs {
#        input:
#            files = [ CallSVs.pbsv_vcf, CallSVs.sniffles_vcf, CallSVs.svim_vcf, CallSVs.cutesv_vcf ],
#            outdir = outdir + "/variants"
#    }

    call FF.FinalizeToDir as FinalizeSmallVariants {
        input:
            files = [ CallSmallVariants.longshot_vcf, CallSmallVariants.longshot_tbi ],
            outdir = outdir + "/variants"
    }

    call FF.FinalizeToDir as FinalizeMergedRuns {
        input:
            files = [ bam, bai ],
            outdir = outdir + "/alignments"
    }

#    call FF.FinalizeToDir as FinalizeMethylation {
#        input:
#            files = [ freq_tsv ],
#            outdir = outdir + "/" + DIR[0] + "/methylation"
#    }
}
