version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/BAMutils.wdl" as BU

import "../../../tasks/Utility/ONTBamShardResetAndDeduplicate.wdl" as CleanOneShard

workflow DeduplicateAndResetONTAlignedBam {
    meta {
        desciption: "Removes duplicate records from an aligned ONT bam, while resetting the alignment information."
    }

    parameter_meta {
        scatter_scheme: "A txt file holding how to scatter the WGS bam. Example (this example size-balance among the shards): ...\nchr5,chr19\nchr6,chrY,chrM\n..."
    }

    input {
        File  aligned_bam
        File? aligned_bai
        File scatter_scheme
    }

    output {
        File result = Merge.res
    }

    # step 1, split input bam by the T2T size-balanced scheme, then for each (don't forget unmapped) shard
    call BU.ShardAlignedBam { input: aligned_bam = aligned_bam, aligned_bai = aligned_bai, scatter_scheme = scatter_scheme }

    # for each shard, step 2
    scatter (shard_bam in ShardAlignedBam.split_bams) {
        call CleanOneShard.Work as DeShard { input: shard_bam = shard_bam }
    }
    call CleanOneShard.Work as DeUmap { input: shard_bam = ShardAlignedBam.unmapped_reads }

    # step 3 gather
    Array[File] fixed_shards = flatten([[DeUmap.clean_bam], DeShard.clean_bam])
    Array[String] basenames_for_merging = flatten([[DeUmap.basename_of_clean_bam], DeShard.basename_of_clean_bam])
    call Utils.ComputeAllowedLocalSSD as SizeIt { input: intended_gb = 50 + ceil(size(fixed_shards, "GB"))}
    call BU.MergeBamsQuerynameSortedWithPicard as Merge { input:
        qns_bams = fixed_shards, base_names = basenames_for_merging,
        out_prefix = basename(aligned_bam, ".bam") + ".DSP-fixed", num_ssds = SizeIt.numb_of_local_ssd
    }
}
