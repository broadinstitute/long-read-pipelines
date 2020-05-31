version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.15/wdl/tasks/PBUtils.wdl" as PB
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.15/wdl/tasks/Utils.wdl" as Utils
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.15/wdl/tasks/AlignReads.wdl" as AR
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.15/wdl/tasks/AlignedMetrics.wdl" as AM
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.15/wdl/tasks/CallSVs.wdl" as SV
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.15/wdl/tasks/Figures.wdl" as FIG
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.15/wdl/tasks/Finalize.wdl" as FF
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.15/wdl/tasks/CallSmallVariants.wdl" as SMV

workflow PBCLRWholeGenomeSingleFlowcell {
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
        String RG = "@RG\\tID:~{ID}\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        call Utils.ShardLongReads { input: unmapped_files = [ subread_bam ], num_reads_per_split = 2000000 }

        scatter (subreads in ShardLongReads.unmapped_shards) {
            call AR.Minimap2 as AlignChunk {
                input:
                    reads      = [ subreads ],
                    ref_fasta  = ref_fasta,
                    RG         = RG,
                    map_preset = "map-pb"
            }
        }

        call AR.MergeBams as MergeChunks { input: bams = AlignChunk.aligned_bam }

        call AM.AlignedMetrics as PerFlowcellSubRunMetrics {
            input:
                aligned_bam    = MergeChunks.merged_bam,
                aligned_bai    = MergeChunks.merged_bai,
                ref_fasta      = ref_fasta,
                ref_dict       = ref_dict,
                ref_flat       = ref_flat,
                dbsnp_vcf      = dbsnp_vcf,
                dbsnp_tbi      = dbsnp_tbi,
                metrics_locus  = metrics_locus,
                per            = "flowcell",
                type           = "subrun",
                label          = ID,
                gcs_output_dir = outdir + "/" + DIR
        }

#        call FIG.Figures as PerFlowcellSubRunFigures {
#            input:
#                summary_files  = [ summary_file ],
#
#                per            = "flowcell",
#                type           = "subrun",
#                label          = SID,
#
#                gcs_output_dir = outdir + "/" + DIR
#        }
    }

    call AR.MergeBams as MergeRuns { input: bams = MergeChunks.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }

    call AM.AlignedMetrics as PerFlowcellRunMetrics {
        input:
            aligned_bam    = MergeRuns.merged_bam,
            aligned_bai    = MergeRuns.merged_bai,
            ref_fasta      = ref_fasta,
            ref_dict       = ref_dict,
            ref_flat       = ref_flat,
            dbsnp_vcf      = dbsnp_vcf,
            dbsnp_tbi      = dbsnp_tbi,
            metrics_locus  = metrics_locus,
            per            = "flowcell",
            type           = "run",
            label          = ID[0],
            gcs_output_dir = outdir + "/" + DIR[0]
    }

#    call FIG.Figures as PerFlowcellRunFigures {
#        input:
#            summary_files  = FindSequencingSummaryFiles.summary_files,
#
#            per            = "flowcell",
#            type           = "run",
#            label          = ID[0],
#
#            gcs_output_dir = outdir + "/" + DIR[0]
#    }

    call SV.CallSVs as CallSVs {
        input:
            bam               = MergeRuns.merged_bam,
            bai               = MergeRuns.merged_bai,

            ref_fasta         = ref_fasta,
            ref_fasta_fai     = ref_fasta_fai,
            tandem_repeat_bed = tandem_repeat_bed
    }

    call SMV.CallSmallVariants as CallSmallVariants {
        input:
            bam               = MergeRuns.merged_bam,
            bai               = MergeRuns.merged_bai,

            ref_fasta         = ref_fasta,
            ref_fasta_fai     = ref_fasta_fai,
            ref_dict          = ref_dict
    }

    ##########
    # Finalize
    ##########

    call FF.FinalizeToDir as FinalizeSVs {
        input:
            files = [ CallSVs.pbsv_vcf,       CallSVs.pbsv_tbi,
                      CallSVs.sniffles_vcf,   CallSVs.sniffles_tbi,
                      CallSVs.svim_vcf,       CallSVs.svim_tbi ],
            outdir = outdir + "/" + DIR[0] + "/variants"
    }

    call FF.FinalizeToDir as FinalizeSmallVariants {
        input:
            files = [ CallSmallVariants.longshot_vcf, CallSmallVariants.longshot_tbi ],
            outdir = outdir + "/" + DIR[0] + "/variants"
    }

    call FF.FinalizeToDir as FinalizeMergedRuns {
        input:
            files = [ MergeRuns.merged_bam, MergeRuns.merged_bai ],
            outdir = outdir + "/" + DIR[0] + "/alignments"
    }
}
