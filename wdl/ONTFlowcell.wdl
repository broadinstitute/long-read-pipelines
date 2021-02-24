version 1.0

import "tasks/ONTUtils.wdl" as ONT
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/CallSVs.wdl" as SV
import "tasks/CallSmallVariants.wdl" as SMV
import "tasks/Figures.wdl" as FIG
import "tasks/Methylation.wdl" as Meth
import "tasks/Finalize.wdl" as FF

workflow ONTFlowcell {
    input {
        File final_summary
        File sequencing_summary
        File ref_map_file

        String participant_name
        Int num_shards = 50

        String gcs_out_root_dir
    }

    parameter_meta {
        final_summary:      "GCS path to '*final_summary*.txt*' file for basecalled fastq files"
        sequencing_summary: "GCS path to '*sequencing_summary*.txt*' file for basecalled fastq files"
        ref_map_file:       "table indicating reference sequence and auxillary file locations"

        participant_name:   "name of the participant from whom these samples were obtained"
        num_shards:         "[default-valued] number of shards into which fastq files should be batched"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTFlowcell/" + participant_name

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

    call Utils.MergeBams as MergeReads { input: bams = AlignReads.aligned_bam, prefix = "~{participant_name}.~{ID}" }

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

#    File bam = MergeReads.merged_bam
#    File bai = MergeReads.merged_bai
#
#    call FF.FinalizeToDir as FinalizeMergedRuns {
#        input:
#            files = [ bam, bai ],
#            outdir = outdir + "/alignments"
#    }

    output {
#        String gcs_outdir = outdir

        File aligned_bam = MergeReads.merged_bam
        File aligned_bai = MergeReads.merged_bai

#        #Float xml_num_records = SummarizeXMLMetadata.xml_num_records
#        #Float xml_total_length = SummarizeXMLMetadata.xml_total_length
#
#        Float num_records = SummarizeSubreadsPBI.results['reads']
#        Float total_length = SummarizeSubreadsPBI.results['bases']
#        Float raw_yield = SummarizeSubreadsPBI.results['yield']
#
#        Float polymerase_mean = SummarizeSubreadsPBI.results['polymerase_mean']
#        Float polymerase_n50 = SummarizeSubreadsPBI.results['polymerase_n50']
#
#        Float subread_mean = SummarizeSubreadsPBI.results['subread_mean']
#        Float subread_n50 = SummarizeSubreadsPBI.results['subread_n50']
#
#        Float ccs_num_records = SummarizeCCSPBI.results['reads']
#        Float ccs_total_length = SummarizeCCSPBI.results['bases']
#        Float ccs_mean_qual = SummarizeCCSPBI.results['mean_qual']
#        Float ccs_yield = SummarizeCCSPBI.results['yield']
#
#        Float ccs_num_records_q20 = SummarizeCCSQ20PBI.results['reads']
#        Float ccs_total_length_q20 = SummarizeCCSQ20PBI.results['bases']
#        Float ccs_mean_qual_q20 = SummarizeCCSQ20PBI.results['mean_qual']
#        Float ccs_yield_q20 = SummarizeCCSQ20PBI.results['yield']
#
#        Float zmws_input = SummarizeCCSReport.zmws_input
#        Float zmws_pass_filters = SummarizeCCSReport.zmws_pass_filters
#        Float zmws_fail_filters = SummarizeCCSReport.zmws_fail_filters
#        Float zmws_pass_filters_pct = SummarizeCCSReport.zmws_pass_filters_pct
#        Float zmws_fail_filters_pct = SummarizeCCSReport.zmws_fail_filters_pct
    }
}
