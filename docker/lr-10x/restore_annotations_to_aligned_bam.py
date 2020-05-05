import argparse
import sys
import time
import os

import tqdm
import pysam

import sqlite3
import base64
import pickle

def main(bam_name, aligned_bam_name, out_bam_name, ignore_tags=None):

    t_start = time.time()

    # Tags to preserve (for now it isn't used):
    tag_names = []
    pysam.set_verbosity(0)  # silence message about the .bai file not being found

    print("Verifying input files exist...")
    files_ok = True
    for f in [bam_name, aligned_bam_name]:
        if not os.path.exists(f):
            print(f"ERROR: Input file does not exist: {f}")
            files_ok = False
    if not files_ok:
        sys.exit(1)
    print("Input files verified.")

    tags_to_ignore = set()
    if ignore_tags:
        print(f"Ignoring the following tags: {' '.join([t for t in ignore_tags])}")
        for t in ignore_tags:
            tags_to_ignore.add(t)
            if t in tag_names:
                tag_names.remove(t)

    # Set up the sqlite database so we can not run out of memory anymore:
    con = sqlite3.connect("read_tags.db")
    cur = con.cursor()
    cur.execute('''CREATE TABLE read_tags (read TEXT PRIMARY KEY, contig TEXT, position INT, tag_tuple_list BLOB)''')
    con.commit()

    print(f"Reading in tags from: {bam_name}")
    with pysam.AlignmentFile(
            bam_name, "rb", check_sq=False, require_index=False
    ) as bam_file, tqdm.tqdm(
        desc="Progress",
        unit=" read",
        file=sys.stderr
    ) as pbar:
        for i, read in enumerate(bam_file):
            tags_to_keep = []
            if len(tag_names) == 0:
                for t in read.get_tags():
                    # Keep the tags we're not ignoring:
                    if t[0] in tags_to_ignore:
                        continue
                    tags_to_keep.append((t[0], t[1]))
            else:
                # Keep only the tags we specified:
                for tag_name in tag_names:
                    try:
                        tags_to_keep.append((tag_name, read.get_tag(tag_name)))
                    except KeyError:
                        pass
            cur.execute(
                f"INSERT INTO read_tags VALUES ("
                f"'{read.query_name}',"
                f"'{read.reference_name}',"
                f"{read.reference_start},"
                f"'{pickle.dumps(tags_to_keep, 0).decode()}')"
            )
            if i % 20000 == 0:
                con.commit()

            pbar.update(1)

    con.commit()
    print("Tag list assembled.")

    print(f"Applying tags to reads in {aligned_bam_name} and writing to {out_bam_name}")
    with pysam.AlignmentFile(
            aligned_bam_name, "rb", check_sq=False, require_index=False
    ) as aligned_bam_file, tqdm.tqdm(
        desc="Progress",
        unit=" read",
        file=sys.stderr
    ) as pbar:

        # Get our header from the input bam file:
        out_bam_header_dict = aligned_bam_file.header.to_dict()

        # Add our program group to it:
        pg_dict = {
            "ID": f"restore-annotations-0.0.1",
            "PN": "restore-annotations",
            "VN": f"0.0.1",
            "DS": 'Reads an unaligned bam with annotations and an aligned bam containing the same reads.  '
                  'Copies tags on the unaligned bam to reads in the aligned file with the same name.',
            "CL": " ".join(sys.argv),
        }
        if "PG" in out_bam_header_dict:
            out_bam_header_dict["PG"].append(pg_dict)
        else:
            out_bam_header_dict["PG"] = [pg_dict]
        out_header = pysam.AlignmentHeader.from_dict(out_bam_header_dict)

        with pysam.AlignmentFile(out_bam_name, "wb", header=out_header) as out_bam_file:
            for read in aligned_bam_file:
                raw = cur.execute(
                    f"SELECT tag_tuple_list FROM read_tags WHERE read=='{read.query_name}'"
                ).fetchone()

                if raw:
                    tag_tuple_list = pickle.loads(raw[0].encode())
                    for tag_name, tag_value in tag_tuple_list:
                        read.set_tag(tag_name, tag_value)

                out_bam_file.write(read)

                pbar.update(1)
    con.close()
    t_end = time.time()
    print("Done")
    print(f"Elapsed time: {t_end - t_start:2.2f}s")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Reads an unaligned bam with annotations and an aligned bam containing the same reads.'
                    'Copies tags on the unaligned bam to reads in the aligned file with the same name.')
    required_named_args = parser.add_argument_group('required named arguments')
    required_named_args.add_argument('-b', '--bam', help='Bam file from which to read annotations.', required=True)
    required_named_args.add_argument('-a', '--aligned-bam', help='Aligned bam to which to apply annotations.',
                                     required=True)
    required_named_args.add_argument('-o', '--out-name',
                                     help='Output bam file name', required=True)

    optional_named_args = parser.add_argument_group('optional named arguments')
    optional_named_args.add_argument('-i', '--ignore-tags', nargs='+',
                                     help='Tags to ignore.  These tags will not be transferred to the resulting file.',
                                     required=False)

    args = parser.parse_args()

    main(args.bam, args.aligned_bam, args.out_name, args.ignore_tags)
