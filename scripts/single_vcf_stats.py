import argparse
from collections import defaultdict
import gzip

import numpy as np
from text_histogram import histogram


p = argparse.ArgumentParser()
p.add_argument("vcf", nargs="+")
args = p.parse_args()

vcf_paths = args.vcf

c = defaultdict(int)
for vcf_path in vcf_paths:
    print("============================")
    alt_allele_sizes = []
    with gzip.open(vcf_path) if vcf_path.endswith(".gz") else open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            ref = fields[3]
            alt_alleles = fields[4]
            for alt in alt_alleles.split(","):
                c['total alleles'] += 1
                if len(ref) < len(alt):
                    c['insertions'] += 1
                    alt_allele_sizes.append(len(alt))
                elif len(ref) > len(alt):
                    c['deletions'] += 1
                elif len(ref) == 1 and len(alt) == 1:
                    c['snps'] += 1

    for label, value in sorted(c.items(), key=lambda x: (x[1], x[0]), reverse=True):
        print("==> %15s : %9d" % (label, value))

    print("log2(alt allele length) for insertions in " + vcf_path)
    histogram(np.log2(alt_allele_sizes), custbuckets=",".join(["%0.1f" % (k/2.0) for k in range(1, 20)]))
