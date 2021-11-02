version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF
import "tasks/Longbow.wdl" as LONGBOW

import "tasks/Structs.wdl"

workflow MasSeqDemultiplex {

    meta {
        description : "This workflow is designed to pre-process data from the MASSeq v2 protocol.  It will first filter PacBio Sequel IIe reads to keep only the high quality reads.  Then it will create several bam files - one for each model - and populate those files with the reads that best fit those models."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String gcs_input_dir
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/MasSeqDemultiplex"

        # Maximum polymerase read length for reads to be included in the output files from the data being processed.
        # Used here so that the CCS reclamation process can still run in the following workflows, but so the data will
        # be small enough to run through Longbow in a reasonable amount of time.
        Int max_read_length = 60000

        Array[String] models = ["mas10", "mas15"]

        String longbow_docker_version = "us.gcr.io/broad-dsp-lrma/lr-longbow:0.4.3"

        String? sample_name
    }

    parameter_meta {
        gcs_input_dir : "Input folder on GCS in which to search for BAM files to process."
        gcs_out_root_dir : "Root output GCS folder in which to place results of this workflow."

        max_read_length : "[optional] Maximum polymerase read length for reads to be included in the output files from the data being processed.  Used here so that the CCS reclamation process can still run in the following workflows, but so the data will be small enough to run through Longbow in a reasonable amount of time.  (default: 60000)."

        models : "[optional] Models to use for demultiplexing.  (default: [\"mas10\', \"mas15\"])."

        sample_name : "[optional] The name of the sample to associate with the data in this workflow."
    }

    # Version of this workflow.
    String VERSION = "0.3"

    # Create an object to disable preemption.  This should only be used for testing.
    RuntimeAttr disable_preemption = object {
        preemptible_tries:  0
    }

    RuntimeAttr disable_preemption_with_longbow = object {
        mem_gb: 32,
        preemptible_tries:  0,
        docker: longbow_docker_version
    }

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_01_WdlExecutionStartTimestamp { input: }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams as t_02_FindBams { input: gcs_input_dir = gcs_input_dir }

    scatter (reads_bam in t_02_FindBams.ccs_bams) {

        call PB.GetRunInfo as t_03_GetRunInfo { input: subread_bam = reads_bam }

        String SM  = select_first([sample_name, t_03_GetRunInfo.run_info["SM"]])
        String PL  = "PACBIO"
        String PU  = t_03_GetRunInfo.run_info["PU"]
        String DT  = t_03_GetRunInfo.run_info["DT"]
        String ID  = PU
        String DS  = t_03_GetRunInfo.run_info["DS"]
        String DIR = SM + "." + ID

        String RG_subreads  = "@RG\\tID:~{ID}.subreads\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
        String RG_consensus = "@RG\\tID:~{ID}.consensus\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
        String RG_array_elements = "@RG\\tID:~{ID}.array_elements\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        # Shard widely so we can go faster:
        File read_pbi = sub(reads_bam, ".bam$", ".bam.pbi")
        call PB.ShardLongReads as t_04_ShardLongReads {
            input:
                unaligned_bam = reads_bam,
                unaligned_pbi = read_pbi,
                prefix = SM + "_shard",
                num_shards = 2000,
                runtime_attr_override = disable_preemption
        }

        scatter (sharded_reads in t_04_ShardLongReads.unmapped_shards) {

            # 0 - Filter out the kinetics tags from PB files:
            call PB.RemoveKineticsTags as t_05_RemoveKineticsTags {
                input:
                    bam = sharded_reads,
                    prefix = SM + "_kinetics_removed"
            }

            # 1 - filter the reads by the maximum length.
            # NOTE: Because of the reclamation process we only limit length here and allow the following workflows to
            #       perform the actual read quality filtering.
            String fbmrq_prefix = basename(t_05_RemoveKineticsTags.bam_file, ".bam")
            call Utils.Bamtools as t_06_FilterByMaxReadLength {
                input:
                    bamfile = sharded_reads,
                    prefix = fbmrq_prefix + "_good_reads",
                    cmd = "filter",
                    args = '-length "<=' + max_read_length + '"',
                    runtime_attr_override = disable_preemption
            }

            # 2 - split the reads by the model:
            String adis_prefix = basename(t_06_FilterByMaxReadLength.bam_out, ".bam")
            call LONGBOW.Demultiplex as t_07_AssignReadsToModels {
                input:
                    bam = t_06_FilterByMaxReadLength.bam_out,
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

        String base_out_dir = outdir + "/" + DIR + "/" + t_01_WdlExecutionStartTimestamp.timestamp_string

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
                    prefix = SM + "_" + model + ".reads",
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
                    outdir = base_out_dir + "/" + SM + "_" + model,
                    runtime_attr_override = disable_preemption
            }

            # Copy over the metadata to our finalized folders so we can just run our workflow on it:
            call PB.CopyMetadataFilesToNewDir as t_12_CopyMetadataForMas10Reads {
                input:
                    input_gs_path = gcs_input_dir,
                    dest_gs_path = base_out_dir + "/" + SM + "_" + model
            }
        }
    }
}
