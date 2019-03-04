import gzip
import sys

MINIMUM_STR_LENGTH = 3
VALID_NUCLEOTIDES = {'A', 'C', 'G', 'T'}

input_vcf = sys.argv[1]

fopen = gzip.open if input_vcf.endswith(".gz") else open

for line in fopen(input_vcf):
    if line.startswith("#"):
        sys.stdout.write(line)
        continue

    fields = line.split("\t")
    ref = fields[3].upper()
    alt = fields[4].upper()

    if len(alt) < 50:
        continue

    line = "\t".join(fields)

    sys.stdout.write(line)
