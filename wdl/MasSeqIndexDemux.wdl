version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/TranscriptAnalysis/Preprocessing_Tasks.wdl" as TX_PRE
import "tasks/Longbow.wdl" as LONGBOW
import "tasks/Finalize.wdl" as FF

workflow MasSeqIndexDemux {

    meta {
        description : "This workflow will split MAS-seq data that is indexed with a 10bp sequence at the 3' end."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String gcs_input_dir
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/MasSeqIndexDemux"

        String mas_seq_model = "mas15"
        Int max_read_length_bases = 60000

        String? sample_name

        # Set up some meta parameters here so we can adjust for when we want things to go VERY fast:
        Int primary_scatter_width = 300
    }

    parameter_meta {
        gcs_input_dir : "Input folder on GCS in which to search for BAM files to process."
        gcs_out_root_dir : "Root output GCS folder in which to place results of this workflow."

        mas_seq_model : "[optional] The name of the model to use to annotate the data in this workflow (default: mas15)."
        max_read_length_bases : "[optional] Maximum number of bases in a read for that read to be processed by this workflow (default: 60000)."

        sample_name : "[optional] The name of the sample to associate with the data in this workflow."

        primary_scatter_width : "[optional] Width to use for the primary scatter operation on this dataset (default: 300)."
    }

    # Version of this workflow.
    String VERSION = "0.1"

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_01_WdlExecutionStartTimestamp { input: }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams as t_02_FindBams { input: gcs_input_dir = gcs_input_dir }

    # Check here if we found ccs bams or subread bams:
    Boolean use_subreads = t_02_FindBams.has_subreads

    # Make sure we have **EXACTLY** one bam file to run on:
    if (length(t_02_FindBams.ccs_bams) != 1) {
        call Utils.FailWithWarning as t_03_FailOnMultiBamFiles { input: warning = "Error: Multiple BAM files found.  Cannot continue!" }
     }

    # Define some runtime attributes for later:
    RuntimeAttr fast_network_attrs_preemptible = object {
        cpu_cores:  4,
        mem_gb:     32,
        disk_type:  "LOCAL",
        preemptible_tries:  1
    }

    # Alias our bam file so we can work with it easier:
    File reads_bam = t_02_FindBams.ccs_bams[0]

    call PB.GetRunInfo as t_04_GetRunInfo { input: subread_bam = reads_bam }

    String SM  = select_first([sample_name, t_04_GetRunInfo.run_info["SM"]])
    String PL  = "PACBIO"
    String PU  = t_04_GetRunInfo.run_info["PU"]
    String DT  = t_04_GetRunInfo.run_info["DT"]
    String ID  = PU
    String DS  = t_04_GetRunInfo.run_info["DS"]
    String DIR = SM + "." + ID

    File read_pbi = sub(reads_bam, ".bam$", ".bam.pbi")
    call PB.ShardLongReads as t_05_ShardLongReads {
        input:
            unaligned_bam = reads_bam,
            unaligned_pbi = read_pbi,
            prefix = SM + "_shard",
            num_shards = primary_scatter_width,
    }

    scatter (sharded_reads in t_05_ShardLongReads.unmapped_shards) {

        # Filter out the kinetics tags from PB files:
        call PB.RemoveKineticsTags as t_06_RemoveKineticsTags {
            input:
                bam = sharded_reads,
                prefix = SM + "_kinetics_removed"
        }

        # Filter out reads that are too long:
        call Utils.Bamtools as t_07_FilterReadsByLength {
            input:
                bamfile = t_06_RemoveKineticsTags.bam_file,
                prefix = SM + "_reads_of_good_length",
                cmd = "filter",
                args = '-length "<=' + max_read_length_bases + '"'
        }

        # Annotate our CCS Corrected reads:
        call LONGBOW.Annotate as t_08_LongbowAnnotateCCSReads {
            input:
                reads = t_07_FilterReadsByLength.bam_out,
                model = mas_seq_model
        }

        call TX_PRE.DemuxMasSeqDataByIndex as t_09_DemuxMasSeqDataByIndex {
            input:
                array_bam = t_08_LongbowAnnotateCCSReads.annotated_bam,
        }
    }

    #################################################
    #   __  __
    #  |  \/  | ___ _ __ __ _  ___
    #  | |\/| |/ _ \ '__/ _` |/ _ \
    #  | |  | |  __/ | | (_| |  __/
    #  |_|  |_|\___|_|  \__, |\___|
    #                   |___/
    #
    # Merge all the sharded files we created above into files for this
    # input bam file.
    #################################################

    call Utils.MergeBams as t_10_MergeDemuxI1 { input: bams = t_09_DemuxMasSeqDataByIndex.demux_i1, prefix = SM + "_demux_i1" }
    call Utils.MergeBams as t_11_MergeDemuxI2 { input: bams = t_09_DemuxMasSeqDataByIndex.demux_i2, prefix = SM + "_demux_i2" }
    call Utils.MergeBams as t_12_MergeDemuxI3 { input: bams = t_09_DemuxMasSeqDataByIndex.demux_i3, prefix = SM + "_demux_i3" }
    call Utils.MergeBams as t_13_MergeDemuxI4 { input: bams = t_09_DemuxMasSeqDataByIndex.demux_i4, prefix = SM + "_demux_i4" }
    call Utils.MergeBams as t_14_MergeDemuxAmbiguous { input: bams = t_09_DemuxMasSeqDataByIndex.demux_ambiguous, prefix = SM + "_demux_ambiguous" }
    call TX_PRE.MergeDemuxMasSeqByIndexLogs as t_15_MergeDemuxMasSeqByIndexLogs { input: demux_logs = t_09_DemuxMasSeqDataByIndex.log_file }

    ######################################################################
    #             _____ _             _ _
    #            |  ___(_)_ __   __ _| (_)_______
    #            | |_  | | '_ \ / _` | | |_  / _ \
    #            |  _| | | | | | (_| | | |/ /  __/
    #            |_|   |_|_| |_|\__,_|_|_/___\___|
    #
    ######################################################################

    # NOTE: We key all finalization steps on the static report.
    #       This will prevent incomplete runs from being placed in the output folders.
    String base_out_dir = outdir + "/" + DIR + "/" + t_01_WdlExecutionStartTimestamp.timestamp_string

        call FF.FinalizeToDir as t_16_FinalizeOutputFiles {
        input:
            files = [
                t_10_MergeDemuxI1.merged_bam,
                t_10_MergeDemuxI1.merged_bai,
                t_11_MergeDemuxI2.merged_bam,
                t_11_MergeDemuxI2.merged_bai,
                t_12_MergeDemuxI3.merged_bam,
                t_12_MergeDemuxI3.merged_bai,
                t_13_MergeDemuxI4.merged_bam,
                t_13_MergeDemuxI4.merged_bai,
                t_14_MergeDemuxAmbiguous.merged_bam,
                t_14_MergeDemuxAmbiguous.merged_bai,
                t_15_MergeDemuxMasSeqByIndexLogs.merged_log,
            ],
            outdir = base_out_dir + "/",
            keyfile = t_15_MergeDemuxMasSeqByIndexLogs.merged_log
    }

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good
    call FF.WriteCompletionFile as t_17_WriteCompletionFile {
        input:
            outdir = base_out_dir + "/",
            keyfile =  t_15_MergeDemuxMasSeqByIndexLogs.merged_log
    }
}
