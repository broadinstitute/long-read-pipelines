version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/BAMutils.wdl" as BU

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

    # task 1, split input bam by the T2T size-balanced scheme, then for each (don't forget unmapped) shard
    call BU.ShardAlignedBam { input: aligned_bam = aligned_bam, aligned_bai = aligned_bai, scatter_scheme = scatter_scheme }

    # for each shard
    scatter (shard_bam in ShardAlignedBam.split_bams) {
    #    task 2, samtools reset and querysort by picard
    #    task 3, deduplicate
        call BU.SamtoolsReset as Magic { input: bam = shard_bam }
        call BU.QuerynameSortBamWithPicard as SortUnaligned { input: bam = Magic.res }
        call BU.DeduplicateQuerynameSortedBam as DeQS { input: qnorder_bam = SortUnaligned.qnsort_bam}
        String base_for_each = basename(DeQS.dedup_bam)
        call BU.GetDuplicateReadnamesInQnameSortedBam as CheckDedupShard { input: qns_bam = DeQS.dedup_bam }
        if ( CheckDedupShard.result_may_be_corrupted || 0!=length(read_lines(CheckDedupShard.dup_names_txt)) ) {
            call Utils.StopWorkflow as DedupShardFail { input: reason = "Deduplication isn't successful for ~{shard_bam}."}
        }
    }
    # and that includes the unmapped reads
    call BU.SamtoolsReset as RemoveDict { input: bam = ShardAlignedBam.unmapped_reads }  # consistent processing with the mapped ones (practically, take away the @SQ lines in the header)
    call BU.QuerynameSortBamWithPicard as QNSortUnammpedReads { input: bam = RemoveDict.res }
    call BU.DeduplicateQuerynameSortedBam as DeU { input: qnorder_bam = QNSortUnammpedReads.qnsort_bam }
    String base_for_unmap = basename(DeU.dedup_bam)
    call BU.GetDuplicateReadnamesInQnameSortedBam as CheckDedupUnmapped { input: qns_bam = DeU.dedup_bam }
    if ( CheckDedupUnmapped.result_may_be_corrupted || 0!=length(read_lines(CheckDedupUnmapped.dup_names_txt)) ) {
        call Utils.StopWorkflow as DedupUnmappedFail { input: reason = "Deduplication isn't successful for ~{DeU.dedup_bam}."}
    }

    # task 4 gather
    Array[File] fixed_shards = flatten([[DeU.dedup_bam], DeQS.dedup_bam]) # takein the unmapped "shard"
    Array[String] basenames_for_merging = flatten([[base_for_unmap], base_for_each])
    call Utils.ComputeAllowedLocalSSD as SizeIt { input: intended_gb = 50 + ceil(size(fixed_shards, "GB"))}
    call BU.MergeBamsQuerynameSortedWithPicard as Merge { input:
        qns_bams = fixed_shards, base_names = basenames_for_merging,
        out_prefix = basename(aligned_bam, ".bam") + ".DSP-fixed", num_ssds = SizeIt.numb_of_local_ssd
    }
}
