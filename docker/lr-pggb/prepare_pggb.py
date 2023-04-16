#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import bz2
import sys
from pathlib import Path

import skbio


def open_compressed(fpath: Path, *args, **kwargs):
    func = {
        ".gz": gzip.open,
        ".bz2": bz2.open,
    }.get(fpath.suffix, open)

    return func(fpath, *args, **kwargs)


def main():
    parser = argparse.ArgumentParser(description='Prepare a single FASTA with multiple genomes to be used as input '
                                                 'for PGGB')

    parser.add_argument(
        'refs', metavar='REF_FASTA', nargs='+', type=Path,
        help="Input reference genomes. Can be compressed with gzip or bzip2."
    )

    parser.add_argument(
        '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
        help="Output FASTA. Defaults to stdout."
    )

    parser.add_argument(
        '-d', '--directory-as-id', action="store_true",
        help="Use the directory name of the FASTA file as reference name instead of the filename."
    )

    parser.add_argument(
        '-D', '--delimiter', default="#", required=False,
        help="Pangenome sequencing naming delimiter. Default: %(default)s. For more information visit "
             "https://github.com/pangenome/PanSN-spec"
    )

    args = parser.parse_args()
    delim = args.delimiter

    for ref in args.refs:
        with open_compressed(ref, "rt") as f:
            if ref.suffix in {".gz", ".bz2"}:
                ref = ref.with_name(ref.stem)

            if args.directory_as_id:
                ref_id = ref.parent.name
            else:
                ref_id = ref.stem

            for r in skbio.io.read(f, "fasta", verify=False):
                # Currently only supports haploid genomes
                r.metadata['id'] = f"{ref_id}{delim}1{delim}{r.metadata['id']}"
                skbio.io.write(r, "fasta", into=args.output)


if __name__ == '__main__':
    main()
