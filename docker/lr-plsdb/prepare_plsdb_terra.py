#!/usr/bin/env python3
import argparse
import functools

print_tsv = functools.partial(print, sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('plsdb_tsv', type=argparse.FileType('r'), help="PLSDB metadata TSV file path")
    parser.add_argument('gcs_plsdb_dir', help="GCS base directory with PLSDB sequences.")

    args = parser.parse_args()

    gcs_base_dir = args.gcs_plsdb_dir
    if not gcs_base_dir.endswith('/'):
        gcs_base_dir += "/"

    # We swap the first two columns of the original TSV (using accession as plasmid ID), and add GCS FASTA URL
    first = True
    for line in args.plsdb_tsv:
        line = line.strip()
        if not line:
            continue

        components = line.split('\t')

        if first:
            print_tsv("entity:plasmid_id", "plasmid_fasta", "UID_NUCCORE", *components[2:])
            first = False
        else:
            gcs_url = f"{gcs_base_dir}{components[1]}/{components[1]}.fna.bz2"
            print_tsv(components[1], gcs_url, components[0], *components[2:])
