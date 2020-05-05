#!/usr/bin/env python3

import os
import sys
import argparse

import pysam
import tqdm


def remove_interval_from_cigar_tuples(cigartuples, interval_start, interval_end):
    """Remove cigar elements from the pysam cigartuples object that correspond
    to the bases in the given interval defined by start and end.

    NOTE: interval_start and interval_end are INCLUSIVE
          and are assumed to be properly ordered (interval_start < interval_end)"""

    out_cigar_tuples = []

    current_element_pos = 0
    for element, length in cigartuples:

        # We have 11 cases for how the interval can relate to the cigar element:
        #                  [   CIGAR   ]
        # 1  -  [     ]                                    - Interval outside Cigar (before)
        # 2  -           [     ]                           - Interval starts before, ends inside Cigar
        # 3  -                 [     ]                     - Interval contained within Cigar
        # 4  -                       [     ]               - Interval starts inside, ends after Cigar
        # 5  -                                  [     ]    - Interval outside Cigar (after)
        # 6  -       [                       ]             - Interval starts before ends after Cigar
        # 7  -             [       ]                       - Interval start equals cigar interval start, end within cigar
        # 8  -             [                ]              - Interval start equals cigar interval start, end outside cigar
        # 9  -          [              ]                   - Interval starts before cigar, end equals cigar interval end
        # 10 -                [        ]                   - Interval starts after cigar interval start, end equals cigar interval end
        # 11 -             [           ]                   - Interval equals cigar interval
        new_element_length = length
        # 1:
        if interval_end < current_element_pos:
            # Cigar doesn't change:
            pass
        # 2:
        elif interval_start < current_element_pos < interval_end < current_element_pos + length:
            new_element_length -= interval_end - current_element_pos
        # 3:
        elif current_element_pos < interval_start < interval_end < current_element_pos + length:
            new_element_length -= interval_end - interval_start
        # 4:
        elif current_element_pos < interval_start < current_element_pos + length < interval_end:
            new_element_length -= current_element_pos + length - interval_start
        # 5:
        elif interval_start > current_element_pos + length:
            # Cigar doesn't change:
            pass
        # 6:
        elif interval_start < current_element_pos < current_element_pos + length < interval_end:
            # Cigar element is completely removed.
            new_element_length = 0
        # 7:
        elif interval_start == current_element_pos and current_element_pos < interval_end < current_element_pos + length:
            new_element_length -= interval_end - current_element_pos - 1
        # 8:
        elif interval_start == current_element_pos and current_element_pos + length < interval_end:
            # Cigar element is completely removed.
            new_element_length = 0
        # 9:
        elif interval_start < current_element_pos and current_element_pos + length == interval_end:
            # Cigar element is completely removed.
            new_element_length = 0
        # 10:
        elif current_element_pos < interval_start < current_element_pos + length and \
                (current_element_pos + length) == interval_end:
            new_element_length -= current_element_pos + length - interval_start - 1
        # 11:
        elif interval_start == current_element_pos and current_element_pos + length == interval_end:
            # Cigar element is completely removed.
            new_element_length = 0
        else:
            # This should never happen:
            raise Exception("This should never happen.")

        # Add in our new cigar element if it still has bases in it:
        if new_element_length > 0:
            out_cigar_tuples.append((element, length))

        # We track the original length because we produce a brand-new cigar based on our original lengths:
        current_element_pos += length

    return out_cigar_tuples


def read_has_good_tags(read):
    for t in ["CB", "XB", "ZU", "XU"]:
        if not read.has_tag(t) or read.get_tag(t) == ".":
            return False
    return True


def main(input_bam, out_bam_name):

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
            "ID": f"extract-cbc-and-umi-from-annotated-read-0.0.1",
            "PN": "extract-cbc-and-umi-from-annotated-read",
            "VN": f"0.0.1",
            "DS": "Extracts the cell barcode and UMI sequences from the given annotaed bam file and appends them to "
                  "the read name.  Base qualities for those reads will be removed as well.  "
                  "The output can be converted into a FASTQ file that is compatible with UMI Tools.",
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

                # Make sure we can actually get the info we need from the read:
                if not read_has_good_tags(read):
                    pbar.update(1)
                    continue

                cbc = read.get_tag("CB")
                cbc_pos = [int(p) for p in read.get_tag("XB").split(',')]
                umi = read.get_tag("ZU")
                umi_pos = [int(p) for p in read.get_tag("XU").split(',')]

                # Append CBC and UMI to read name:
                read_name = f"{read.query_name}_{cbc}_{umi}"

                # Get our quals and bases:
                bases = read.query_sequence
                quals = read.query_qualities
                # cigartuples = read.cigartuples

                # Remove the quals and bases of our UMI and CBC:
                # NOTE: We remove from the back to the front of the read
                # so we don't have to worry about changing indices.
                if cbc_pos[0] > umi_pos[0]:
                    bases = bases[0:cbc_pos[0]] + bases[cbc_pos[1]:]
                    bases = bases[0:umi_pos[0]] + bases[umi_pos[1]:]

                    quals = quals[0:cbc_pos[0]] + quals[cbc_pos[1]:]
                    quals = quals[0:umi_pos[0]] + quals[umi_pos[1]:]

                    # cigartuples = remove_interval_from_cigar_tuples(cigartuples, cbc_pos[0], cbc_pos[1])
                    # cigartuples = remove_interval_from_cigar_tuples(cigartuples, umi_pos[0], umi_pos[1])
                else:
                    bases = bases[0:umi_pos[0]] + bases[umi_pos[1]:]
                    bases = bases[0:cbc_pos[0]] + bases[cbc_pos[1]:]

                    quals = quals[0:umi_pos[0]] + quals[umi_pos[1]:]
                    quals = quals[0:cbc_pos[0]] + quals[cbc_pos[1]:]

                    # cigartuples = remove_interval_from_cigar_tuples(cigartuples, umi_pos[0], umi_pos[1])
                    # cigartuples = remove_interval_from_cigar_tuples(cigartuples, cbc_pos[0], cbc_pos[1])

                # Create an output record to write to the bam file:
                read.query_name = read_name
                read.query_sequence = bases
                read.query_qualities = quals

                # Remove the tags that give position information of the CBC and UMI:
                # NOTE: We do this because these bases are excised from the reads.
                read.set_tag("XB", None)
                read.set_tag("XU", None)

                # read.cigartuples = cigartuples

                # Write out our read:
                out_bam_file.write(read)

                pbar.update(1)
                reads_processed += 1
    print("Done!")
    print(f"Reads processed: {reads_processed}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Extracts the cell barcode and UMI sequences from the given annotaed bam file and appends them to "
                    "the read name.  Base qualities for those reads will be removed as well.  "
                    "The output can be converted into a FASTQ file that is compatible with UMI Tools.",
        epilog="The input file must have been annotated with the 10x tool.  "
               "Each read must contain the CB, XB, ZU, and XU tags (as defined in tool.py).  "
               "NOTE: Assumes that the given bam file is not aligned."
    )

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-b', '--bam',
                               help='Unaligned bam file from which to extract the cell barcode and UMI sequences.',
                               required=True)
    requiredNamed.add_argument('-o', '--out-name',
                               help='Output bam file name',
                               required=True)

    args = parser.parse_args()
    main(args.bam, args.out_name)
