version 1.0

import "Utils.wdl"
import "BAMutils.wdl" as BU

workflow Work {

    meta {
        desciption:
        "Given a shard of an aligned ONT bam (sharded by chromosome), reset alignment and remove duplicate."
        note:
        "This is purely an artifact because WDL files cannot have two workflows."
    }
    input {
        File shard_bam
    }
    output {
        File clean_bam = DeQS.dedup_bam
        String basename_of_clean_bam = basename(DeQS.dedup_bam)
    }

    call BU.SamtoolsReset as Magic { input: bam = shard_bam }
    call BU.QuerynameSortBamWithPicard as SortUnaligned { input: bam = Magic.res }
    call BU.DeduplicateQuerynameSortedBam as DeQS { input: qnorder_bam = SortUnaligned.qnsort_bam}

    # verify
    call BU.GetDuplicateReadnamesInQnameSortedBam as InitialCheckDedupShard { input: qns_bam = DeQS.dedup_bam, trial_idx = 1 }
    if ( InitialCheckDedupShard.result_may_be_corrupted ) {
        call BU.GetDuplicateReadnamesInQnameSortedBam as RetryCheckDedupShard { input: qns_bam = DeQS.dedup_bam, trial_idx = 2 }
        if ( RetryCheckDedupShard.result_may_be_corrupted ) {
            call BU.GetDuplicateReadnamesInQnameSortedBam as LastCheckDedupShard { input: qns_bam = DeQS.dedup_bam, trial_idx = 3, trial_max = 3 }
        }
    }

    # do not change order
    Boolean CheckOperationFailed = select_first([LastCheckDedupShard.result_may_be_corrupted,
                                                 RetryCheckDedupShard.result_may_be_corrupted,
                                                 InitialCheckDedupShard.result_may_be_corrupted])
    Array[String] dup_names = read_lines(select_first([LastCheckDedupShard.dup_names_txt,
                                                       RetryCheckDedupShard.dup_names_txt,
                                                       InitialCheckDedupShard.dup_names_txt]))
    if ( CheckOperationFailed || 0!=length(dup_names) ) {
        call Utils.StopWorkflow as DedupShardFail { input: reason = "Deduplication isn't successful for ~{shard_bam}."}
    }
}
