#!/usr/bin/env python3
"""Generate a synthetic VCF exercising every hmmibd-rs --from-bcf filter.

Sites are crafted so each filter branch is hit at least once, plus a 2-cluster
clonal structure on the 'good' sites so IBD inference produces real segments.
GT is deliberately set to ./. everywhere to prove read_dom ignores GT and
re-derives the call from FORMAT/AD.
"""
NSAMP = 12
samples = [f"S{i:02d}" for i in range(NSAMP)]
GROUP0 = range(0, 6)   # clonal cluster A
GROUP1 = range(6, 12)  # clonal cluster B

HDR = [
    "##fileformat=VCFv4.2",
    '##FILTER=<ID=PASS,Description="pass">',
    '##FILTER=<ID=LowQual,Description="low quality">',
    '##contig=<ID=1,length=1000000>',
    '##contig=<ID=2,length=1000000>',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="allele depths">',
]
cols = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"] + samples


def ad_for(nallele, dom):
    """depth vector: 20 on the dominant allele, 0 elsewhere (clear call)."""
    v = [0]*nallele
    v[dom] = 20
    return v


def ad_missing(nallele):
    """sub-threshold total depth -> read_dom emits -1 (missing)."""
    v = [0]*nallele
    v[0] = 2
    if nallele > 1:
        v[1] = 1
    return v


def fmt_cell(adv):
    return "./.:" + ",".join(str(x) for x in adv)


def row(chrom, pos, ref, alt, filt, per_sample_ad):
    cells = [fmt_cell(a) for a in per_sample_ad]
    return "\t".join([chrom, str(pos), ".", ref, alt, ".", filt, ".", "GT:AD"] + cells)


rows = []

# ---- chr1: 25 good biallelic SNP sites, 2-cluster clonal structure ----
for k in range(25):
    pos = 500*(k+1)
    a0 = k % 2            # cluster A's allele this site
    a1 = 1 - a0           # cluster B's allele
    ad = []
    for i in range(NSAMP):
        allele = a0 if i in GROUP0 else a1
        # inject a few deterministic singletons/missing for informative variety
        if (k*7 + i) % 23 == 0:
            ad.append(ad_missing(2))            # occasional missing
        elif (k*5 + i) % 17 == 0:
            ad.append(ad_for(2, 1 - allele))    # occasional flip
        else:
            ad.append(ad_for(2, allele))
    rows.append((1, pos, row("1", pos, "A", "C", "PASS", ad)))

# close-pair: 2bp from the pos=500 site -> dropped by min_snp_sep (default 5)
rows.append((1, 502, row("1", 502, "A", "C", "PASS",
             [ad_for(2, 0 if i in GROUP0 else 1) for i in range(NSAMP)])))

# multiallelic SNP (3 SNP alleles), PASS -> KEPT (maf = 2nd-allele / total)
ad = []
for i in range(NSAMP):
    allele = i // 4          # 0,0,0,0,1,1,1,1,2,2,2,2
    ad.append(ad_for(3, allele))
rows.append((1, 13000, row("1", 13000, "A", "C,G", "PASS", ad)))

# non-PASS site -> dropped (filter_column_pass_only)
rows.append((1, 13500, row("1", 13500, "A", "C", "LowQual",
             [ad_for(2, 0 if i in GROUP0 else 1) for i in range(NSAMP)])))

# indel major/minor -> dropped (major_minor_alleles_must_be_snps); ALT=AT is indel
rows.append((1, 14000, row("1", 14000, "A", "AT", "PASS",
             [ad_for(2, 0 if i in GROUP0 else 1) for i in range(NSAMP)])))

# monomorphic (all allele 0) -> maf 0 < min_maf -> dropped
rows.append((1, 14500, row("1", 14500, "A", "C", "PASS",
             [ad_for(2, 0) for i in range(NSAMP)])))

# low-depth everywhere -> all missing -> non_missing_rate 0 < min_site_nonmissing
rows.append((1, 15000, row("1", 15000, "A", "C", "PASS",
             [ad_missing(2) for i in range(NSAMP)])))

# ---- chr2: 15 good biallelic SNP sites ----
for k in range(15):
    pos = 1000*(k+1)
    a0 = (k+1) % 2
    a1 = 1 - a0
    ad = []
    for i in range(NSAMP):
        allele = a0 if i in GROUP0 else a1
        if (k*3 + i) % 19 == 0:
            ad.append(ad_missing(2))
        else:
            ad.append(ad_for(2, allele))
    rows.append((2, pos, row("2", pos, "A", "C", "PASS", ad)))

# coordinate-sort (chrom as int, then pos) and emit
rows.sort(key=lambda r: (r[0], r[1]))
with open("test.vcf", "w") as fh:
    fh.write("\n".join(HDR) + "\n")
    fh.write("\t".join(cols) + "\n")
    fh.write("\n".join(r[2] for r in rows) + "\n")
print(f"wrote {len(rows)} records, {NSAMP} samples")
