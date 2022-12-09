version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF
import "tasks/Longbow.wdl" as LONGBOW

workflow LongbowCorrect {

    meta {
        description : "This workflow runs Longbow Correct on a given bam file."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File reads_bam
        String sample_name

        String mas_seq_model = "mas_15+sc_10x5p"

        Int ccs_lev_dist = 2
        Int clr_lev_dist = 3

        String raw_tag = "CR"
        String final_tag = "CB"

        File cell_barcode_whitelist = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/737K-august-2016.txt"

        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/LongbowCorrect"
    }

    parameter_meta {
        reads_bam : "Input file to process."
        sample_name : "Name of the sample in the given reads file."
        mas_seq_model : "[optional] built-in mas-seq model to use (Default: mas_15+sc_10x5p)"
        gcs_out_root_dir : "[optional] Root output GCS folder in which to place results of this workflow.  (Default: gs://broad-dsde-methods-long-reads-outgoing/LongbowCorrect)"
    }

    ##############################################################################################################

    call Utils.GetCurrentTimestampString as t_001_WdlExecutionStartTimestamp { input: }

    call PB.PBIndex as t_002_PbIndexReads {
        input:
            bam = reads_bam
    }

    call PB.ShardLongReads as t_003_ShardLongReads {
        input:
            unaligned_bam = reads_bam,
            unaligned_pbi = t_002_PbIndexReads.pbindex,
            prefix = sample_name + "_shard",
            num_shards = 100,
    }

    scatter (main_shard_index in range(length(t_003_ShardLongReads.unmapped_shards))) {
        File sharded_reads = t_003_ShardLongReads.unmapped_shards[main_shard_index]
        call LONGBOW.Correct as t_004_CorrectReads {
            input:
                reads = sharded_reads,
                barcode_allow_list = cell_barcode_whitelist,
                model = mas_seq_model,
                ccs_lev_dist_threshold = ccs_lev_dist,
                clr_lev_dist_threshold = clr_lev_dist,
                prefix = sample_name + "_longbow_corrected_shard_" + main_shard_index,
                raw_barcode_tag = raw_tag,
                corrected_barcode_tag = final_tag,
        }
    }

    call Utils.MergeBams as t_005_MergeCorrectedReads {
        input:
            bams = t_004_CorrectReads.corrected_barcodes_bam,
            prefix = sample_name + "_longbow_corrected"
    }

    call Utils.MergeBams as t_006_MergeUnCorrectedReads {
        input:
            bams = t_004_CorrectReads.uncorrected_barcodes_bam,
            prefix = sample_name + "_longbow_uncorrected"
    }


    call PB.PBIndex as t_007_PbIndexCorrectedReads {
        input:
            bam = t_005_MergeCorrectedReads.merged_bam
    }

    call PB.PBIndex as t_008_PbIndexUncorrectedReads {
        input:
            bam = t_006_MergeUnCorrectedReads.merged_bam
    }


    ##############################################################################################################

    String outdir = sub(gcs_out_root_dir, "/$", "")
    String base_out_dir = outdir + "/" + sample_name + "/" + t_001_WdlExecutionStartTimestamp.timestamp_string

    # Finalize the final corrected, aligned array elements:
    call FF.FinalizeToDir as t_009_FinalizeQuantifiedArrayElements {
        input:
            files = [
                t_005_MergeCorrectedReads.merged_bam,
                t_005_MergeCorrectedReads.merged_bai,
                t_007_PbIndexCorrectedReads.pbindex,

                t_006_MergeUnCorrectedReads.merged_bam,
                t_006_MergeUnCorrectedReads.merged_bai,
                t_008_PbIndexUncorrectedReads.pbindex,
            ],
            outdir = base_out_dir
    }

}