version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF
import "tasks/Longbow.wdl" as LONGBOW

workflow LongbowAnnotate {

    meta {
        description : "This workflow runs Longbow Annotate on a given bam file."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File reads_bam
        String sample_name

        String mas_seq_model = "mas15"

        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/LongbowAnnotate"
    }

    parameter_meta {
        reads_bam : "Input file to process."
        sample_name : "Name of the sample in the given reads file."
        mas_seq_model : "[optional] built-in mas-seq model to use (Default: mas15)"
        gcs_out_root_dir : "[optional] Root output GCS folder in which to place results of this workflow.  (Default: gs://broad-dsde-methods-long-reads-outgoing/LongbowAnnotate)"
    }

    ##############################################################################################################

    call Utils.GetCurrentTimestampString as t_01_WdlExecutionStartTimestamp { input: }

    call PB.PBIndex as t_02_PbIndexReads {
        input:
            bam = reads_bam
    }

    call PB.ShardLongReads as t_03_ShardLongReads {
        input:
            unaligned_bam = reads_bam,
            unaligned_pbi = t_02_PbIndexReads.pbindex,
            prefix = sample_name + "_shard",
            num_shards = 300,
    }

    RuntimeAttr new_longbow_docker = object {
        docker: "us.gcr.io/broad-dsp-lrma/lr-longbow:0.3.1"
    }

    scatter (sharded_reads in t_03_ShardLongReads.unmapped_shards) {
        call LONGBOW.Annotate as t_04_AnnotateReads {
            input:
                reads = sharded_reads,
                model = mas_seq_model,
                runtime_attr_override = new_longbow_docker
        }
    }

    call Utils.MergeBams as t_05_MergeAnnotatedReads {
        input:
            bams = t_04_AnnotateReads.annotated_bam,
            prefix = sample_name + "_longbow_annotated"
    }

    call PB.PBIndex as t_06_PbIndexReads {
        input:
            bam = t_05_MergeAnnotatedReads.merged_bam
    }


    ##############################################################################################################

    String outdir = sub(gcs_out_root_dir, "/$", "")
    String base_out_dir = outdir + "/" + sample_name + "/" + t_01_WdlExecutionStartTimestamp.timestamp_string

    # Finalize the final annotated, aligned array elements:
    call FF.FinalizeToDir as t_07_FinalizeQuantifiedArrayElements {
        input:
            files = [
                t_05_MergeAnnotatedReads.merged_bam,
                t_05_MergeAnnotatedReads.merged_bai,
                t_06_PbIndexReads.pbindex,
            ],
            outdir = base_out_dir
    }

}