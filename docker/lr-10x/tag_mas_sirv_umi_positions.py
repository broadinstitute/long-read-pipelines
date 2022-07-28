#!/usr/bin/env python3

import os
import sys
import argparse

import pysam
import tqdm

gDEFAULT_UMI_LENGTH_bases = 8
gUMI_TAG_NAME = "ZU"
gANNOTATED_SEQUENCES_TAG = "SG"


def main(input_bam, out_bam_name, umi_length, umi_tag_name):

    global gANNOTATED_SEQUENCES_TAG

    print("Verifying input file(s) exist...", file=sys.stderr)
    files_ok = True
    for f in [input_bam]:
        if not os.path.exists(f):
            print(f"ERROR: Input file does not exist: {f}", file=sys.stderr)
            files_ok = False
    if not files_ok:
        sys.exit(1)
    print("Input files verified.", file=sys.stderr)

    reads_processed = 0
    reads_filtered = 0

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
            "ID": f"tag-mas-sirv-umi-positions-0.0.1",
            "PN": "tag-mas-sirv-umi-positions",
            "VN": f"0.0.1",
            "DS": "Annotates the position of the UMI in each read.",
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

                # Get the tag for the 10x adapter:
                sequences = read.get_tag(gANNOTATED_SEQUENCES_TAG).split(",")
                found_10x_adapter = False
                for s in sequences:
                    # Get the tag name and the coordinates:
                    name, coords = s.split(":")
                    if name == "10x_Adapter":
                        found_10x_adapter = True
                        coords = [int(c) for c in coords.split('-')]
                        break

                # Get the UMI coords:
                if found_10x_adapter:
                    umi_start = coords[1] + 1
                    umi_end = umi_start + umi_length - 1  # Subtract 1 because of inclusive coordinates.
                else:
                    # Read does not have a 10x adapter.  If it starts with a "random" block and the read name specifies
                    # that the first segment is "START" then we assume the 10x adapter was not present and pull the
                    # bases from the first random block.
                    if sequences[0].startswith("random:") and "/START-" in read.query_name:
                        umi_start = 0
                        umi_end = umi_start + umi_length - 1
                    else:
                        # There is no UMI in this read:
                        print(f"WARNING: Read does not contain 10x adapter - cannot find UMI: {read.query_name}", file=sys.stderr)
                        reads_filtered += 1
                        continue

                # Get the UMI sequence from the read:
                umi_seq = read.query_sequence[umi_start:umi_end+1]

                # Set our tag for our UMI:
                read.set_tag(umi_tag_name, umi_seq)

                # Adjust our random tag:
                for i, s in enumerate(sequences):
                    if s.startswith("random:"):
                        coords = [int(c) for c in s.split(':')[1].split("-")]
                        # Adjust start position:
                        coords[0] += umi_length
                        break
                sequences[i] = f"random:{coords[0]}-{coords[1]}"
                read.set_tag(gANNOTATED_SEQUENCES_TAG, ",".join(sequences))

                # Adjust our sequence and bases:
                bases = read.query_sequence[0:umi_start] + read.query_sequence[umi_end+1:]
                quals = read.query_qualities[0:umi_start] + read.query_qualities[umi_end + 1:]

                read.query_sequence = bases
                read.query_qualities = quals

                # Write out our read:
                out_bam_file.write(read)

                # Update our counts:
                pbar.update(1)
                reads_processed += 1
    print("Done!", file=sys.stderr)
    print(f"Reads processed: {reads_processed}", file=sys.stderr)
    print(f"Reads filtered:  {reads_filtered}", file=sys.stderr)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=f"Extracts the UMI from each read in the given bam file.  "
                    f"UMI Position starts at the base immediately after the end of the 10x-adapter (TSO) "
                    f"and is {gDEFAULT_UMI_LENGTH_bases} bases long.  "
                    f"NOTE: Assumes that the given file contains only annotated MAS-seq SIRV reads.",
        epilog=f"The UMI is added as a tag ({gUMI_TAG_NAME}) to each read.  The sequence data and the 'random' tag "
               f"containing the UMI sequence itself are adjusted to reflect the bases assigned to the UMI."
    )

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-b', '--bam',
                               help='Aligned bam file from which to extract the contig information into a tag.',
                               required=True)
    requiredNamed.add_argument('-o', '--out-name',
                               help='Output bam file name',
                               required=True)

    args = parser.parse_args()
    main(args.bam, args.out_name, gDEFAULT_UMI_LENGTH_bases, gUMI_TAG_NAME)
