#!/usr/bin/env python3
import bz2
import argparse
from pathlib import Path

import skbio


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('plsdb_fasta', type=Path,
                        help="Path to PLSDB fasta file")
    parser.add_argument('output_dir', type=Path, default=None,
                        help="Path to output directory. Defaults to same directory as the FASTA file.")

    args = parser.parse_args()

    output_dir = args.output_dir if args.output_dir is not None else args.plsdb_fasta.parent
    total = 0
    with bz2.open(args.plsdb_fasta, "rt") as f:
        for r in skbio.io.read(f, "fasta"):
            plasmid_id = r.metadata['id']
            plasmid_dir = output_dir / plasmid_id
            plasmid_dir.mkdir(exist_ok=True, parents=True)
            print("Writing", plasmid_id)
            with bz2.open(output_dir / plasmid_id / f"{plasmid_id}.fna.bz2", "wt") as o:
                skbio.io.write(r, "fasta", into=o)

            total += 1

    print("Written", total, "plasmid sequences.")
