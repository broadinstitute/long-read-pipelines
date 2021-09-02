version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF
import "tasks/AlignReads.wdl" as AR

workflow Minimap2 {

    meta {
        description : "Align the given reads with Minimap2"
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String input_bam

        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/Minimap2"

        File ref_fasta =  "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
        File ref_fasta_index = "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai"
        File ref_fasta_dict = "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict"

        String minimap2_preset = "asm20"

        String? sample_name
    }

    # Create some runtime attributes that will force google to do network transfers really fast:
    RuntimeAttr fast_network_attrs = object {
        cpu_cores:  4,
        mem_gb:     32,
        disk_type:  "LOCAL",
        preemptible_tries:  0
    }
    RuntimeAttr fast_network_attrs_preemptible = object {
        cpu_cores:  4,
        mem_gb:     32,
        disk_type:  "LOCAL",
        preemptible_tries:  1
    }

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_01_WdlExecutionStartTimestamp { input: }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    File read_pbi = sub(input_bam, ".bam$", ".bam.pbi")
    call PB.ShardLongReads as t_02_ShardLongReads {
        input:
            unaligned_bam = input_bam,
            unaligned_pbi = read_pbi,
            prefix = sample_name + "_shard",
            num_shards = 300,
    }

    scatter (sharded_reads in t_02_ShardLongReads.unmapped_shards) {

        ## No more preemption on this sharding - takes too long otherwise.
        RuntimeAttr disable_preemption_runtime_attrs = object {
            preemptible_tries: 0
        }

        call AR.Minimap2 as t_03_AlignReads {
            input:
                reads      = [ sharded_reads ],
                ref_fasta  = ref_fasta,
                map_preset = minimap2_preset
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

    call Utils.MergeBams as t_04_MergeAlignedReads { input: bams = t_03_AlignReads.aligned_bam, prefix = sample_name + "_aligned_reads" }
    call PB.PBIndex as t_05_PbIndexReads {
        input:
            bam = t_04_MergeAlignedReads.merged_bam
    }

    ######################################################################
    #             _____ _             _ _
    #            |  ___(_)_ __   __ _| (_)_______
    #            | |_  | | '_ \ / _` | | |_  / _ \
    #            |  _| | | | | | (_| | | |/ /  __/
    #            |_|   |_|_| |_|\__,_|_|_/___\___|
    #
    ######################################################################

    String base_out_dir = outdir + "/" + sample_name + "/" + t_01_WdlExecutionStartTimestamp.timestamp_string

    ##############################################################################################################
    # Finalize the final annotated, aligned array elements:
    call FF.FinalizeToDir as t_06_FinalizeAlignedReads {
        input:
            files = [
                t_04_MergeAlignedReads.merged_bam,
                t_04_MergeAlignedReads.merged_bai,
                t_05_PbIndexReads.pbindex,
            ],
            outdir = base_out_dir + "/",
    }

    ##############################################################################################################
    # Write out completion file so in the future we can be 100% sure that this run was good:
    call FF.WriteCompletionFile as t_07_WriteCompletionFile {
        input:
            outdir = base_out_dir + "/",
            keyfile = t_05_PbIndexReads.pbindex
    }
}
