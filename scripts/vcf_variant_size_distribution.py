import gzip
import os
import sys

import numpy as np
from text_histogram import histogram

input_vcf = sys.argv[1]

fopen = gzip.open if input_vcf.endswith(".gz") else open

alt_allele_sizes = []
for line in fopen(input_vcf):
    if line.startswith("#"):
        continue
    
    fields = line.split("\t")
    #ref = fields[3]
    len_alt = len(fields[4])
    
    if len_alt > 1:
        alt_allele_sizes.append(len_alt)

print("-----")
print("log2(alt allele length) for " + input_vcf)
histogram(np.log2(alt_allele_sizes), custbuckets=",".join(["%0.1f" % (k/2.0) for k in range(1, 20)]))
