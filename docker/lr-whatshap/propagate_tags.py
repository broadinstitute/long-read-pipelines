#!/usr/bin/env python

import argparse
import pysam


def main():
    parser = argparse.ArgumentParser(description='Propagate BAM tags', prog='propagate_tags')
    parser.add_argument('-t', '--tags', type=str, nargs='+', action='append', help="BAM tags")
    parser.add_argument('-s', '--source', type=str, help='source BAM')
    parser.add_argument('-d', '--dest', type=str, help='destination BAM')
    parser.add_argument('-o', '--output', type=str, help='output BAM')
    args = parser.parse_args()

    srtags = dict()

    print("Storing tags...")
    pysam.set_verbosity(0)
    with pysam.AlignmentFile(args.source, "rb", check_header=False, check_sq=False) as s:
        for sr in s:
            tags = dict()
            for tag in args.tags:
                if sr.has_tag(tag[0]):
                    tags[tag[0]] = sr.get_tag(tag[0])

            srtags[sr.query_name] = tags

    print(f"Stored tags for {len(srtags)} reads.")

    print("Applying tags...")
    num_tagged = 0
    with pysam.AlignmentFile(args.dest, "rb", check_header=False, check_sq=False) as d:
        with pysam.Samfile(args.output, 'wb', header=d.header) as out:
            for dr in d:
                is_tagged = False

                for tag in args.tags:
                    if dr.query_name in srtags and tag[0] in srtags[dr.query_name]:
                        old_tag = dr.get_tag(tag[0]) if dr.has_tag(tag[0]) else None
                        dr.set_tag(tag[0], srtags[dr.query_name][tag[0]], replace=True)
                        is_tagged = True

                        print(f'{dr.query_name} {tag[0]} {old_tag} {dr.get_tag(tag[0])}')

                if is_tagged:
                    num_tagged += 1

                out.write(dr)

    print(f"Tagged {num_tagged} reads.")


if __name__ == "__main__":
    main()
