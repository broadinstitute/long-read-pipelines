version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF
import "tasks/Longbow.wdl" as LONGBOW

import "tasks/Structs.wdl"

workflow MasSeqDemultiplexTriage {

    meta {
        description : "This workflow is designed to pre-process data from the MASSeq v2 protocol.  It will first filter PacBio Sequel IIe reads to keep only the high quality reads.  Then it will create several bam files - one for each model - and populate those files with the reads that best fit those models."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        Array[File] input_shards
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/MasSeqDemultiplex"

        Array[String] models = ["mas10", "mas15"]

        String longbow_docker_version = "us.gcr.io/broad-dsp-lrma/lr-longbow:0.4.3"

        String sample_name
        Int chunk
    }

    # Version of this workflow.
    String VERSION = "0.3"

    # Create an object to disable preemption.  This should only be used for testing.
    RuntimeAttr disable_preemption = object {
        preemptible_tries:  0
    }

    RuntimeAttr disable_preemption_with_longbow = object {
        preemptible_tries:  0,
        docker: longbow_docker_version
    }

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_01_WdlExecutionStartTimestamp { input: }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    scatter (reads_bam in input_shards) {

        # 2 - split the reads by the model:
        String adis_prefix = basename(reads_bam, ".bam") + "_chunk-" + chunk
        call LONGBOW.Demultiplex as t_07_AssignReadsToModels {
            input:
                bam = reads_bam,
                prefix = adis_prefix,
                models = models,
                runtime_attr_override = disable_preemption_with_longbow
        }
    }

    # Create an object to make things faster by allocating more resources:
    RuntimeAttr bigger_resources_for_network = object {
        cpu_cores:  4,
        mem_gb:     8,
        preemptible_tries:  0
    }

    String base_out_dir = outdir + "/TRIAGE/" + sample_name + "/" + t_01_WdlExecutionStartTimestamp.timestamp_string

    # NOTE: We have to do something stupid here because WDL can't seem to iterate through the nested array
    #       the way that I need it to
    Array[File] all_demuxed_bams = flatten(t_07_AssignReadsToModels.demultiplexed_bams)

    # Scatter by each model so we can consolidate files from each one:
    scatter ( model in models ) {

        call Utils.FilterListOfStrings as t_08_FilterDemuxedBamsByModel {
            input :
                list_to_filter = all_demuxed_bams,
                query = "_" + model + ".bam$"
            }

        call Utils.MergeBams as t_09_MergeModelBams {
            input:
                bams = t_08_FilterDemuxedBamsByModel.filtered_list,
                prefix = sample_name + "_" + model + ".reads",
                runtime_attr_override = bigger_resources_for_network
        }
        call PB.PBIndex as t_10_PbIndexModelBam {
            input:
                bam = t_09_MergeModelBams.merged_bam,
                runtime_attr_override = bigger_resources_for_network
        }

        call FF.FinalizeToDir as t_11_FinalizeMergedModelBam {
            input:
                files = [
                    t_09_MergeModelBams.merged_bam,
                    t_09_MergeModelBams.merged_bai,
                    t_10_PbIndexModelBam.pbindex,
                ],
                outdir = base_out_dir + "/" + sample_name + "_" + model,
                runtime_attr_override = disable_preemption
        }
    }
}
