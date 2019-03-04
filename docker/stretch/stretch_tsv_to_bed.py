"""
Parses the STRetch output .tsv format (described in https://github.com/Oshlack/STRetch/wiki/Output-files):

chrom
start
end
sample
repeatunit
reflen - number of repeat units in the reference
locuscoverage - number of STR reads assigned to that locus
outlier - z score testing for outliers
p_adj - adjusted p value, is this locus significantly expanded relative to other samples? (p values have already been adjusted for multiple testing using the Benjamini-Hochberg method)
bpInsertion - estimated size of allele in bp inserted relative to the reference
repeatUnits - estimated total size of allele in repeat units
"""

import argparse

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-o", "---output-bed-path")
    p.add_argument("stretch_tsv_path")
    args = p.parse_args()

    if not args.stretch_tsv_path.endswith(".tsv"):
        p.error("Invalid input filename: %s" % args.stretch_tsv_path)

    return args


def main():

    args = parse_args()

    if not args.output_bed_path:
        args.output_bed_path = args.stretch_tsv_path.replace(".tsv", "") + ".bed"

    with open(args.stretch_tsv_path) as input_file, open(args.output_bed_path, "w") as output_file:
        for i, line in enumerate(input_file):
            line = line.rstrip("\n")
            fields = line.split("\t")
            if i == 0:
                # header should be:
                # chrom   start   end     sample  repeatunit      reflen  locuscoverage   outlier p_adj   bpInsertion     repeatUnits
                header = fields
                assert header[0] == "chrom"
                continue

            record = dict(zip(header, fields))

            output_file.write("\t".join(map(str, [
                record["chrom"],
                int(record["start"]) - 1,  # convert to 0-based coordinates in .bed
                int(record["end"]),
                record["repeatunit"],
                float(record["repeatUnits"]) - float(record["reflen"]),  # how many repeat units are in the variant = total repeat units - repeat units in reference
                ";".join([
                    f"p_adj={record['p_adj']}",
                    f"z_score={record['outlier']}",
                    f"coverage={record['locuscoverage']}"
                ]),
            ])) + "\n")

    print(f"Wrote {i} lines to {args.output_bed_path}")

if __name__ == "__main__":
    main()
