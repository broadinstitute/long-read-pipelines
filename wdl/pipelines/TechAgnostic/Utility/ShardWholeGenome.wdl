version 1.0

import "../../../tasks/Utility/Utils.wdl"

import "../../../tasks/Utility/BAMutils.wdl" as BU

workflow Split {
    meta {
        description: "Split input BAM aligned to a reference genome."
    }
    parameter_meta {
        contig_filter: "List of contigs in the ref genomes to skip. Ignored when ref_scatter_interval_list* are provided (assuming users know how they want to shard)."
        ref_scatter_interval_list_locator: "A file holding paths to interval_list files; provide when explicit sharding scheme is desired."
        ref_scatter_interval_list_ids: "A file that gives short IDs to the interval_list files; provide when explicit sharding scheme is desired."

        id_bam_bai_of_shards: "Id, bam and bai (in the form of <id, <bam, bai>>) of the sharded intput BAM"
    }
    input {
        File ref_dict
        File bam
        File bai
        Array[String] contig_filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']
        File? ref_scatter_interval_list_locator
        File? ref_scatter_interval_list_ids
    }

    output {
        Array[Pair[String, Pair[File, File]]] id_bam_bai_of_shards = zip(ids_of_interval_lists, sharded_bam_bais)
    }

    if (defined(ref_scatter_interval_list_locator)) { # custom
        File scatter_interval_list_ids = select_first([ref_scatter_interval_list_ids])
        File scatter_interval_list_loc = select_first([ref_scatter_interval_list_locator])
        Array[String] custom_interval_list_ids   = read_lines(scatter_interval_list_ids)
        Array[String] custom_interval_list_files = read_lines(scatter_interval_list_loc)
        Array[Pair[String, String]] ided_interval_list_files = zip(custom_interval_list_ids, custom_interval_list_files)

        scatter (pair in ided_interval_list_files) {

            call Utils.ResilientSubsetBam as user_controled_split {
                input:
                    bam = bam,
                    bai = bai,
                    interval_list_file = pair.right,
                    interval_id = pair.left,
                    prefix = basename(bam, ".bam")
            }
            if (user_controled_split.is_samtools_failed) {
                # attempt again, first retry streaming,
                call Utils.ResilientSubsetBam as retry_streaming {
                    input:
                        bam = bam,
                        bai = bai,
                        interval_list_file = pair.right,
                        interval_id = pair.left,
                        prefix = basename(bam, ".bam")
                }
                # then localize if that still fails
                if (retry_streaming.is_samtools_failed) {
                    call BU.SubsetBamToLocusLocal {
                        input:
                            bam = bam,
                            bai = bai,
                            interval_list_file = pair.right,
                            interval_id = pair.left,
                            prefix = basename(bam, ".bam")
                    }
                }
                # call Utils.StopWorkflow as CustomShardFailedStreaming { input: reason = "Streaming from BAM for subsetting into ~{pair.left} failed."}
            }
            File custom_shard_bam = select_first([SubsetBamToLocusLocal.subset_bam, retry_streaming.subset_bam, user_controled_split.subset_bam])
            File custom_shard_bai = select_first([SubsetBamToLocusLocal.subset_bai, retry_streaming.subset_bai, user_controled_split.subset_bai])
        }
    }

    if (!defined(ref_scatter_interval_list_locator)) { # per contig/chromosome
        call Utils.MakeChrIntervalList {
            input:
                ref_dict = ref_dict,
                filter = contig_filter
        }
        scatter (c in MakeChrIntervalList.chrs) {
            String contig = c[0]
            call Utils.SubsetBam as default_split {
                input:
                    bam = bam,
                    bai = bai,
                    locus = contig
            }
            if (default_split.is_samtools_failed) {
                call Utils.StopWorkflow as StandardShardFailedStreaming { input: reason = "Streaming from BAM for subsetting into ~{contig} failed."}
            }
            String default_interval_list_id = contig
        }
    }

    Array[Pair[File, File]] sharded_bam_bais = zip(select_first([custom_shard_bam, default_split.subset_bam]),
                                                   select_first([custom_shard_bai, default_split.subset_bai])
                                                  )
    Array[String] ids_of_interval_lists = select_first([default_interval_list_id, custom_interval_list_ids])
}
