import subprocess
import re
import sys

from google.cloud import storage
import xml.etree.ElementTree as etree
import datetime
import argparse
import pprint
import os

parser = argparse.ArgumentParser(description='Detect long read run information.')
parser.add_argument('--SM', type=str, help='sample name (SM) to use in place of detected value')
parser.add_argument('--ID', type=str, help='read group identifier (ID) to use in place of detected value')
parser.add_argument('--CN', type=str, help='center name (CN) to use in place of the detected value')
parser.add_argument('--PL', type=str, help='platform name (PL) to use in place of the detected value')
parser.add_argument('--PM', type=str, help='platform model (PM) to use in place of the detected value')
parser.add_argument('--PU', type=str, help='platform unit (PU) to use in place of the detected value')
parser.add_argument('--DS', type=str, help='description (DS) to use in place of the detected value')
parser.add_argument('gcs_path', metavar='S', type=str, help='GCS path for long read data to process')
args = parser.parse_args()

delimiter = '/'
gcs_path = args.gcs_path.replace("gs://", "").split(delimiter)

bucket_name = gcs_path[0]
prefix = delimiter.join(gcs_path[1:])

info = {}
bams = []
fqs = []
f5s = []

ri = {
    'ID': None,
    'CN': 'BI',
    'DT': datetime.datetime.now().isoformat(),
    'PL': None,
    'PM': None,
    'PU': None,
    'SM': None,
    'DS': None,
    'SO': "NA",
    'DA': None,
}

if prefix.endswith(".bam") and not prefix.endswith(".scraps.bam"):
    bams.append(prefix)
elif prefix.endswith(".fastq.gz"):
    fqs.append(prefix)
else:
    if prefix[-1] != delimiter:
        prefix += delimiter

    storage_client = storage.Client()
    blobs = storage_client.list_blobs(bucket_name, prefix=prefix)

    pp = pprint.PrettyPrinter(indent=4)

    # First, populate run info (ri) hash with information gleaned from run metadata
    for blob in blobs:
        if blob.name.endswith(".bam") and not blob.name.endswith(".scraps.bam"):
            bams.append(blob.name)

        if not "fail" in blob.name and bool(re.search('(f(ast)?q)(.gz)?', blob.name)):
            fqs.append(blob.name)

        #if not "fail" in blob.name and bool(re.search('(f(ast)?5)(.gz)?', blob.name)):
        if not "fail" in blob.name and blob.name.endswith(".fast5"):
            f5s.append(blob.name)

        if blob.name.endswith("final_summary.txt"):
            blob.download_to_filename("final_summary.txt")

            with open("final_summary.txt") as fp:
                for line in fp:
                    f = line.rstrip().split("=")
                    info[f[0]] = f[1]

            ri['DT'] = info['started']
            ri['PL'] = 'ONT'
            ri['PM'] = info['instrument']
            ri['PU'] = info['flow_cell_id'] + "." + info['position']
            ri['SM'] = info['sample_id']
            ri['ID'] = info['flow_cell_id'] + "." + info['position'] + "." + info['sample_id']
            ri['SO'] = 'unsorted'

        if blob.name.endswith(".metadata.xml"):
            blob.download_to_filename("metadata.xml")

            a = etree.parse("metadata.xml")
            for e in a.iter():
                if bool(e.attrib):
                    tag = re.sub(r"{.+}", "", e.tag)
                    for k, v in e.attrib.items():
                        key = tag + "." + k
                        info[key] = "".join(v)

            movieName = info['CollectionMetadata.Context']
            readType = "SUBREAD"

            ri['ID'] = f'{movieName}//{readType}'
            ri['DT'] = info['WellSample.CreatedAt']
            ri['PL'] = 'PACBIO'
            ri['PM'] = info['CollectionMetadata.InstrumentName']
            ri['PU'] = info['CollectionMetadata.Context']
            ri['SM'] = 'unknown'

            if 'BioSample.Name' in info:
                ri['SM'] = info['BioSample.Name']
            elif 'WellSample.Name' in info:
                ri['SM'] = info['WellSample.Name']

if len(bams) > 0:
    for name in bams:
        bam = f'gs://{bucket_name}/{name}'
        o = subprocess.check_output(["samtools", "view", "-H", bam])
        os = str(o).split("\\n")

        for line in list(filter(lambda x: '@RG' in x, os)):
            for piece in line.split("\\t"):
                if ":" in piece:
                    (k, v) = piece.split(":", maxsplit=1)
                    ri[k] = v

        for line in list(filter(lambda x: '@HD' in x, os)):
            if "coordinate" in line:
                ri['SO'] = "coordinate"
            elif "queryname" in line:
                ri['SO'] = "queryname"
            elif "unsorted" in line:
                ri['SO'] = "unsorted"
            elif "unknown" in line:
                ri['SO'] = "unknown"
elif len(fqs) > 0:
    name = fqs[0]
    fq = f'gs://{bucket_name}/{name}'

    o = ""

    if fq.endswith(".gz"):
        o = subprocess.check_output(f'gsutil cat {fq} | gunzip -c | head -1', shell=True)
    else:
        o = subprocess.check_output(f'gsutil cat {fq} | head -1', shell=True)

    if not 'PL' in ri:
        if '/' in str(o):
            ri['PL'] = 'PACBIO'
        else:
            ri['PL'] = 'ONT'
    if not 'ID' in ri:
        ri['ID'] = os.path.basename(re.sub('.(f(ast)?q)(.gz)?', '', name))
    if not 'SO' in ri:
        ri['SO'] = "unsorted"

if len(bams) > 0:
    ri['DA'] = ",".join(map(lambda x: f'gs://{bucket_name}/{x}', bams))
    ri['TY'] = "BAM"
elif len(fqs) > 0:
    ri['DA'] = ",".join(map(lambda x: f'gs://{bucket_name}/{x}', fqs))
    ri['TY'] = "FASTQ"

if len(fqs) > 0:
    ri['F5'] = ",".join(map(lambda x: f'gs://{bucket_name}/{x}', f5s))

if args.SM is not None:
    ri['SM'] = args.SM
if args.ID is not None:
    ri['ID'] = args.ID
if args.CN is not None:
    ri['CN'] = args.CN
if args.PL is not None:
    ri['PL'] = args.PL
if args.PM is not None:
    ri['PM'] = args.PM
if args.PU is not None:
    ri['PU'] = args.PU
if args.DS is not None:
    ri['DS'] = args.DS

for key in ri.keys():
    ri[key] = ri[key][0] if isinstance(ri[key], tuple) else ri[key]

for k, v in ri.items():
    print(f'{k}\t{v}')

if 0 < len(f5s) != len(fqs):
    sys.exit(f"Error: possible incomplete upload; len(fast5s) ({len(f5s)} != len(fastqs) ({len(fqs)})")
if ri['SM'] is None:
    sys.exit("Error: read group SM not auto-detected nor provided.")
if ri['ID'] is None:
    sys.exit("Error: read group ID not auto-detected nor provided.")
