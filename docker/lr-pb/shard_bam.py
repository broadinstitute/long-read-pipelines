import argparse
import gzip
from collections import OrderedDict
from math import ceil
import pysam
from construct import *

from multiprocessing.pool import ThreadPool
from functools import partial


def write_shard(bam, sharding_offsets, zmw_counts_exp, prefix, nocheck, index):
    """
    Write subset of PacBio bam to a shard, taking care not to split reads
    from the same ZMW across separate files.  These shards are thus suitable
    for correction via CCS.
    """

    bf = pysam.Samfile(bam, 'rb', check_sq=False)

    # Advance to the specified virtual file offset.
    bf.seek(sharding_offsets[index])

    num_reads = 0
    zmw_counts_act = {}
    with pysam.Samfile(f'{prefix}{index}.bam', 'wb', header=bf.header) as out:
        # Write until we've advanced to (but haven't written) the read that begins the next shard.
        while True:
            read = bf.__next__()
            out.write(read)

            # Count the ZMW numbers seen.
            zmw = read.get_tag("zm")
            zmw_counts_act[zmw] = zmw_counts_act.get(zmw, 0) + 1
            num_reads += 1

            if bf.tell() >= sharding_offsets[index+1]:
                break

    # Verify that the count for ZMWs written to this shard match the ZMW count determined
    # from the original read of the index file.  If this exception is thrown, it may indicate
    # that reads from the same ZMW have been erroneously sharded to separate files.
    if not nocheck:
        if zmw_counts_act[zmw] != zmw_counts_exp[zmw]:
            raise Exception(f'Number of reads from a specific ZMW mismatches between the original data '
                            f'and the sharded data ({zmw}: {zmw_counts_exp[zmw]} != {zmw_counts_act[zmw]})')

    return num_reads


def compute_shard_offsets(pbi_file, num_shards):
    """
    Compute all possible shard offsets (keeping adjacent reads from the ZMW together)
    """

    # Decode PacBio .pbi file.  This is not a full decode of the index, only the parts we need for sharding.
    # More on index format at https://pacbiofileformats.readthedocs.io/en/9.0/PacBioBamIndex.html .

    fmt = Struct(
        # Header
        "magic" / Const(b"PBI\x01"),
        "version_patch" / Int8ul,
        "version_minor" / Int8ul,
        "version_major" / Int8ul,
        "version_empty" / Int8ul,
        "pbi_flags" / Int16ul,
        "n_reads" / Int32ul,
        "reserved" / Padding(18),

        # Basic information section (columnar format)
        "rgId" / Padding(this.n_reads * 4),
        "qStart" / Padding(this.n_reads * 4),
        "qEnd" / Padding(this.n_reads * 4),
        "holeNumber" / Array(this.n_reads, Int32sl),
        "readQual" / Padding(this.n_reads * 4),
        "ctxtFlag" / Padding(this.n_reads * 1),
        "fileOffset" / Array(this.n_reads, Int64sl),
    )

    # Make a list of bgzf virtual file offsets for sharding and store ZMW counts.
    file_offsets_hash = OrderedDict()
    last_offset = 0
    zmw_count_hash = {}
    with gzip.open(pbi_file, "rb") as f:
        idx_contents = fmt.parse_stream(f)

        for j in range(0, idx_contents.n_reads):
            # Save only the virtual file offset for the first ZMW hole number, so
            # that shard boundaries always keep reads from the same ZMW together.
            if idx_contents.holeNumber[j] not in file_offsets_hash:
                file_offsets_hash[idx_contents.holeNumber[j]] = idx_contents.fileOffset[j]

            last_offset = idx_contents.fileOffset[j]

            zmw_count_hash[idx_contents.holeNumber[j]] = zmw_count_hash.get(idx_contents.holeNumber[j], 0) + 1

    file_offsets = list(file_offsets_hash.values())
    shard_offsets = []
    for j in range(0, len(file_offsets), ceil(len(file_offsets) / num_shards)):
        shard_offsets.append(file_offsets[j])

    # For the last read in the file, pad the offset so the final comparison in write_shard() retains the final read.
    offset_padding = 100
    shard_offsets.append(last_offset + offset_padding)

    return shard_offsets, zmw_count_hash, idx_contents.n_reads


def main():
    parser = argparse.ArgumentParser(description='Shard .bam file using the .pbi index', prog='shard_bam')
    parser.add_argument('-p', '--prefix', type=str, default="shard", help="Shard filename prefix")
    parser.add_argument('-n', '--num_shards', type=int, default=4, help="Number of shards")
    parser.add_argument('-t', '--num_threads', type=int, default=2, help="Number of threads to use during sharding")
    parser.add_argument('-i', '--index', type=str, required=False, help="PBI index filename")
    parser.add_argument('--nocheck', required=False, help="Do not check for zmw counts.", action="store_true")
    parser.add_argument('bam', type=str, help="BAM")
    args = parser.parse_args()

    pbi = args.bam + ".pbi" if args.index is None else args.index

    if args.nocheck:
        print(f"Disabling ZMW read count integrity check.")

    # Decode PacBio .pbi file and determine the shard offsets.
    print(f"Reading index ({pbi}). This may take a few minutes...", flush=True)
    offsets, zmw_counts, read_count = compute_shard_offsets(pbi, args.num_shards)

    # Prepare a function with arguments partially filled in.
    func = partial(write_shard, args.bam, offsets, zmw_counts, args.prefix, args.nocheck)
    idx = list(range(0, len(offsets) - 1))

    # Write the shards using the specified number of threads.
    print(f"Writing {len(idx)} shards using {args.num_threads} threads...", flush=True)
    res = ThreadPool(args.num_threads).imap_unordered(func, idx)

    # Emit final stats on the sharding.
    all_num_reads_written = list(res)
    count = 0
    for i in range(len(all_num_reads_written)):
        count += all_num_reads_written[i]
        print(f'  - wrote {all_num_reads_written[i]} reads to {args.prefix}{i}.bam', flush=True)

    print(f'Sharded {count}/{read_count} reads across {len(idx)} shards.', flush=True)


if __name__ == "__main__":
    main()
