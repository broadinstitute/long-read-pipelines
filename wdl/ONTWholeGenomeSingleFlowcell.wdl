version 1.0

import "tasks/ONTUtils.wdl" as ONT
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/CallSVs.wdl" as SV
import "tasks/Figures.wdl" as FIG
import "tasks/Finalize.wdl" as FF
import "tasks/CallSmallVariants.wdl" as SMV

workflow ONTWholeGenomeSingleFlowcell {
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

    call ONT.FindSequencingSummaryFiles { input: gcs_input_dir = gcs_input_dir }

    scatter (summary_file in FindSequencingSummaryFiles.summary_files) {
        call ONT.GetRunInfo { input: summary_file = summary_file }
        call ONT.ListFiles as ListFast5s { input: summary_file = summary_file, suffix = "fast5" }
        call ONT.ListFiles as ListFastqs { input: summary_file = summary_file, suffix = "fastq" }

        String SM  = select_first([sample_name, GetRunInfo.run_info["sample_id"]])
        String PL  = "ONT"
        String PU  = GetRunInfo.run_info["instrument"]
        String DT  = GetRunInfo.run_info["started"]
        String ID  = GetRunInfo.run_info["flow_cell_id"] + "." + GetRunInfo.run_info["position"]
        String DIR = GetRunInfo.run_info["protocol_group_id"] + "." + SM + "." + ID
        String SID = ID + "." + sub(GetRunInfo.run_info["protocol_run_id"], "-.*", "")
        String RG = "@RG\\tID:~{SID}\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        call ONT.PartitionManifest as PartitionFast5Manifest { input: manifest = ListFast5s.manifest, N = 4  }
        call ONT.PartitionManifest as PartitionFastqManifest { input: manifest = ListFastqs.manifest, N = 50 }

        scatter (manifest_chunk in PartitionFastqManifest.manifest_chunks) {
            call AR.Minimap2 as AlignChunk {
                input:
                    reads      = read_lines(manifest_chunk),
                    ref_fasta  = ref_fasta,
                    RG         = RG,
                    map_preset = "map-ont"
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
                label          = SID,
                gcs_output_dir = outdir + "/" + DIR
        }

        call FIG.Figures as PerFlowcellSubRunFigures {
            input:
                summary_files  = [ summary_file ],

                per            = "flowcell",
                type           = "subrun",
                label          = SID,

                gcs_output_dir = outdir + "/" + DIR
        }
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

    call FIG.Figures as PerFlowcellRunFigures {
        input:
            summary_files  = FindSequencingSummaryFiles.summary_files,

            per            = "flowcell",
            type           = "run",
            label          = ID[0],

            gcs_output_dir = outdir + "/" + DIR[0]
    }

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
