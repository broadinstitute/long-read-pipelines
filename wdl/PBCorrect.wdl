version 1.0

# Copyright Broad Institute, 2019

import "CorrectReads.wdl" as CR
import "Utils.wdl" as Utils
import "Finalize.wdl" as FF

workflow LRWholeGenomeSingleSample {
    input {
        Array[String] gcs_dirs
        String sample_name
        String gcs_output_dir
    }

    String outdir = sub(sub(gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")

    scatter (gcs_dir in gcs_dirs) {
        call Utils.DetectRunInfo as DetectRunInfo {
            input:
                gcs_dir = gcs_dir,
                sample_name = sample_name,
        }

        call Utils.PrepareRun as PrepareRun { input: files = DetectRunInfo.files, }

        call SLR.ShardLongReads as ShardLongReads { input: unmapped_bam = PrepareRun.unmapped_bam, }

        scatter (unmapped_shard in ShardLongReads.unmapped_shards) {
            call CR.CCS as CCS {
                input:
                    unmapped_shard = unmapped_shard,
                    platform = DetectRunInfo.run_info['PL'],
            }
        }

        call MB.MergeBams as MergeCorrected {
            input:
                aligned_shards = CCS.shard,
                merged_name="~{sample_name}.corrected.unaligned.bam",
        }
    }

    call MB.MergeBams as MergeAllCorrected {
        input:
            aligned_shards = MergeCorrected.merged,
            merged_name="~{sample_name}.corrected.unaligned.bam",
    }

    call VB.ValidateBam as ValidateAllCorrected {
        input:
            input_bam = MergeAllCorrected.merged,
    }

    call FF.FinalizeToDir as FinalizeCorrectedBams {
        input:
            files = [ MergeAllCorrected.merged, MergeAllCorrected.merged_bai ],
            outdir = outdir + "/reads"
    }
}
