import subprocess
import argparse
import pprint
from collections import OrderedDict
import re

parser = argparse.ArgumentParser(description='Detect long read run information.')
parser.add_argument('--SM', type=str, help='sample name (SM) to use in place of detected value')
parser.add_argument('--ID', type=str, help='read group identifier (ID) to use in place of detected value')
parser.add_argument('--CN', type=str, help='center name (CN) to use in place of the detected value')
parser.add_argument('--PL', type=str, help='platform name (PL) to use in place of the detected value')
parser.add_argument('--PM', type=str, help='platform model (PM) to use in place of the detected value')
parser.add_argument('--PU', type=str, help='platform unit (PU) to use in place of the detected value')
parser.add_argument('--DS', type=str, help='description (DS) to use in place of the detected value')
parser.add_argument('unaligned_bam', metavar='U', type=str, help='Unaligned BAM')
args = parser.parse_args()

unaligned_header = subprocess.check_output(f'samtools view -H {args.unaligned_bam} | grep "^@RG"', shell=True)

h = OrderedDict()

uh = str(unaligned_header).split("\\t")
for u in uh:
    if ":" in u:
        k, v = u.split(":", maxsplit=1)
        h[k] = re.sub("[\"']|\\\\n", "", v)

if args.SM is not None:
    h['SM'] = args.SM
if args.ID is not None:
    h['ID'] = args.ID
if args.CN is not None:
    h['CN'] = args.CN
if args.PL is not None:
    h['PL'] = args.PL
if args.PM is not None:
    h['PM'] = args.PM
if args.PU is not None:
    h['PU'] = args.PU
if args.DS is not None:
    h['DS'] = args.DS

newrg = ["@RG"]
for k,v in h.items():
    newrg.append(f'{k}:{v}')

print("\\t".join(newrg), end='')
