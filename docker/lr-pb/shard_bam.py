import argparse
import gzip
from collections import OrderedDict
from math import floor
import pysam
from construct import *

from multiprocessing.pool import ThreadPool
from functools import partial


def write_shard(bam, offsets, prefix, index):
    """
    Write subsection of bam to a shard.
    """

    bf = pysam.Samfile(bam, 'rb', check_sq=False)

    # Advance to the specified virtual file offset.
    bf.seek(offsets[index])

    num_reads = 0
    with pysam.Samfile(f'{prefix}{index}.bam', 'wb', header=bf.header) as out:
        # Write until we've advanced to (but haven't written) the read that begins the next shard.
        while True:
            read = bf.__next__()
            out.write(read)

            num_reads += 1
            if bf.tell() >= offsets[index+1]:
                break

    return num_reads


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Shard .bam file using the .pbi index', prog='shard_bam')
    parser.add_argument('-p', '--prefix', type=str, default="shard", help="Shard filename prefix")
    parser.add_argument('-n', '--num_shards', type=int, default=4, help="Number of shards")
    parser.add_argument('-t', '--num_threads', type=int, default=2, help="Number of threads to use during sharding")
    parser.add_argument('-i', '--index', type=str, required=False, help="Index filename")
    parser.add_argument('bam', type=str, help="BAM")
    args = parser.parse_args()

    pbi = args.bam + ".pbi" if args.index is None else args.index

    print(f"Reading index ({pbi}). This may take a few minutes...")

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
        "rgId" / Array(this.n_reads, Int32sl),
        "qStart" / Array(this.n_reads, Int32sl),
        "qEnd" / Array(this.n_reads, Int32sl),
        "holeNumber" / Array(this.n_reads, Int32sl),
        "readQual" / Array(this.n_reads, Float32b),
        "ctxtFlag" / Array(this.n_reads, Int8ul),
        "fileOffset" / Array(this.n_reads, Int64sl),
    )

    # Make a list of bgzf virtual file offsets for sharding.
    file_offsets_hash = OrderedDict()

    last_offset = 0
    with gzip.open(pbi, "rb") as f:
        a = fmt.parse_stream(f)

        for i in range(0, a.n_reads):
            # Save only the virtual file offset for the first ZMW hole number, so that shard boundaries
            # always keep reads from the same ZMW together.
            if a.holeNumber[i] not in file_offsets_hash:
                file_offsets_hash[a.holeNumber[i]] = a.fileOffset[i]

            last_offset = a.fileOffset[i]

    file_offsets = list(file_offsets_hash.values())

    shard_offsets = []
    for i in range(0, len(file_offsets), floor(len(file_offsets)/args.num_shards)):
        shard_offsets.append(file_offsets[i])

    # For the last read in the file, pad the offset so the final comparison in write_shard() retains the final read.
    shard_offsets.append(last_offset + 100)

    # Prepare a function with arguments partially filled in (the later imap_unordered() call requires functions that
    # only have one remaining argument to be specified).
    func = partial(write_shard, args.bam, shard_offsets, args.prefix)
    idx = list(range(0, len(shard_offsets) - 1))

    # Write the shards using the specified number of threads.
    print(f"Writing {len(idx)} shards using {args.num_threads} threads...")
    res = ThreadPool(args.num_threads).imap_unordered(func, idx)

    # Emit final stats on the sharding.
    all_num_reads_written = list(res)
    count = 0
    for i in range(len(all_num_reads_written)):
        count += all_num_reads_written[i]
        print(f'  - wrote {all_num_reads_written[i]} reads to {args.prefix}{i}.bam')

    assert count == a.n_reads

    print(f'Expected {a.n_reads} reads; sharded {count} reads to {len(idx)} shards.')

