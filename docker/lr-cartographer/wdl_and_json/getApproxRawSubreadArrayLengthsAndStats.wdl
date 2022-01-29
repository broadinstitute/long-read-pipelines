version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Cartographer.wdl" as CART

workflow GetApproxRawSubreadArrayLengthsAndStats {

    input {
        File subreads
        File segments_fasta

        String sample_name

        String? alignment_method = "BWA_ALN"
    }

    ## No more preemption on this sharding - takes too long otherwise.
    RuntimeAttr no_preempt_runtime_attrs = object {
        preemptible_tries: 0
    }

    call Utils.ShardLongReadsWithCopy {
        input:
            unmapped_files = [ subreads ],
            num_reads_per_split = 48000
#            num_reads_per_split = 200000
    }

    scatter (shard in ShardLongReadsWithCopy.unmapped_shards) {

        # Get ZMW Subread stats here to shard them out wider and make it faster:
        call PB.CollectZmwSubreadStats {
            input:
                subreads = shard,
                prefix = sample_name + "_zmw_subread_stats",
                runtime_attr_override = no_preempt_runtime_attrs
        }

        # Get approximate subread array lengths here:
        call CART.GetApproxRawSubreadArrayLengths {
            input:
                reads_file = shard,
                delimiters_fasta = segments_fasta,
                min_qual = 7.0,
                ignore_seqs = ["Poly_A", "Poly_T", "3_prime_TSO", "5_prime_TSO", "TSO"],
                prefix = sample_name + "_approx_raw_subread_array_lengths",
                alignment_algorithm = alignment_method,
                preemptible_attempts = 0
        }
    }

    # Merge our shards of subread stats:
    call Utils.MergeTsvFiles as MergeShardedZmwSubreadStats {
        input:
            tsv_files = CollectZmwSubreadStats.zmw_subread_stats,
            prefix = sample_name + "_zmw_subread_stats"
    }

    # Merge our shards of raw subread array element counts:
    call Utils.MergeTsvFiles as MergeShardedRawSubreadArrayElementCounts {
        input:
            tsv_files = GetApproxRawSubreadArrayLengths.approx_subread_array_lengths,
            prefix = sample_name + "_approx_subread_array_lengths"
    }


    output {
        File zmw_subread_stats            = MergeShardedZmwSubreadStats.merged_tsv
        File approx_subread_array_lengths = MergeShardedRawSubreadArrayElementCounts.merged_tsv
    }
}
