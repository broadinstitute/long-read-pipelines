#!/usr/bin/env python3

import os
import sys
import argparse

import pysam
import tqdm

gGENE_TAG_NAME = "XG"


def main(input_bam, out_bam_name, gene_tag_name):

    print("Verifying input files exist...")
    files_ok = True
    for f in [input_bam]:
        if not os.path.exists(f):
            print(f"ERROR: Input file does not exist: {f}")
            files_ok = False
    if not files_ok:
        sys.exit(1)
    print("Input files verified.")

    reads_processed = 0

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
            input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file, tqdm.tqdm(
        desc="Progress",
        unit=" read",
        file=sys.stderr
    ) as pbar:

        # Get our header from the input bam file:
        out_bam_header_dict = bam_file.header.to_dict()

        # Add our program group to it:
        pg_dict = {
            "ID": f"copy-contig-name-to-tag-0.0.1",
            "PN": "copy-contig-name-to-tag",
            "VN": f"0.0.1",
            "DS": "Copies the name of the contig to which a read is aligned to a tag within that read.",
            "CL": " ".join(sys.argv),
        }
        if "PG" in out_bam_header_dict:
            out_bam_header_dict["PG"].append(pg_dict)
        else:
            out_bam_header_dict["PG"] = [pg_dict]
        out_header = pysam.AlignmentHeader.from_dict(out_bam_header_dict)

        # Open our output file so we can write to it:
        with pysam.AlignmentFile(out_bam_name, "wb", header=out_header) as out_bam_file:
            for read in bam_file:

                # Set our tag for our contig:
                read.set_tag(gene_tag_name, read.reference_name)

                # Write out our read:
                out_bam_file.write(read)

                pbar.update(1)
                reads_processed += 1
    print("Done!")
    print(f"Reads processed: {reads_processed}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=f"Copies the name of the contig to which a read is aligned to a tag within that read.  "
                    f"(tag name: {gGENE_TAG_NAME})",
    )

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-b', '--bam',
                               help='Aligned bam file from which to extract the contig information into a tag.',
                               required=True)
    requiredNamed.add_argument('-o', '--out-name',
                               help='Output bam file name',
                               required=True)

    args = parser.parse_args()
    main(args.bam, args.out_name, gGENE_TAG_NAME)
