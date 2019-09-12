from __future__ import print_function
import argparse
from google.cloud import storage
import xml.etree.ElementTree as etree
import re
import sys


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


parser = argparse.ArgumentParser(description='Create an input JSON for the LRWGS pipeline')
parser.add_argument('--SM', type=str, help='sample name (SM) to use')
parser.add_argument('gcs_paths', metavar='G', type=str, nargs='+', help='GCS path(s)')
args = parser.parse_args()

delimiter = '/'
quoted_paths = []

sample_names = set()
for gcs_path_full in args.gcs_paths:
    quoted_paths.append(f'"{gcs_path_full}"')

    gcs_path = gcs_path_full.replace("gs://", "").split(delimiter)
    bucket_name = gcs_path[0]
    prefix = delimiter.join(gcs_path[1:])

    storage_client = storage.Client()
    blobs = storage_client.list_blobs(bucket_name, prefix=prefix)

    num_seq_blobs = 0
    for blob in blobs:
        if (blob.name.endswith(".bam") and not blob.name.endswith(".scraps.bam")) or bool(re.search('(f(ast)?q)(.gz)?', blob.name)):
            num_seq_blobs += 1

        if blob.name.endswith(".metadata.xml") and not blob.name.endswith(".run.metadata.xml"):
            blob.download_to_filename("metadata.xml")
            a = etree.parse("metadata.xml")
            info = {}

            for e in a.iter():
                if bool(e.attrib):
                    tag = re.sub(r"{.+}", "", e.tag)
                    for k, v in e.attrib.items():
                        key = tag + "." + k
                        info[key] = "".join(v)

            sample_name = info['BioSample.Name'],
            sample_names.add(sample_name[0])

    if num_seq_blobs == 0:
        raise ValueError(f'No sequence data at the location {gcs_path_full}')

eprint(f'[INFO] given= {args.SM} ; detected= {", ".join(sample_names)}')
all_paths = ",\n\t".join(quoted_paths)

s = f"""
    "LRWholeGenomeSingleSample.gcs_dirs": [
        {all_paths}
    ],

    "LRWholeGenomeSingleSample.sample_name": "{args.SM}",

    "LRWholeGenomeSingleSample.ref_fasta":     "gs://broad-dsde-methods-long-reads/resources/references/grch38/Homo_sapiens_assembly38.fasta",
    "LRWholeGenomeSingleSample.ref_fasta_fai": "gs://broad-dsde-methods-long-reads/resources/references/grch38/Homo_sapiens_assembly38.fasta.fai",
    "LRWholeGenomeSingleSample.ref_dict":      "gs://broad-dsde-methods-long-reads/resources/references/grch38/Homo_sapiens_assembly38.dict",

    "LRWholeGenomeSingleSample.tandem_repeat_bed": "gs://broad-dsde-methods-long-reads/resources/references/grch38/human_GRCh38_no_alt_analysis_set.trf.bed",

    "LRWholeGenomeSingleSample.mt_chr_name":   "chrM"
"""

print("{" + s + "}")
