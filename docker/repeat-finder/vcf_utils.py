import vcf

def parse_vcf_file(vcf_path, limit=None):

    reader = vcf.VCFReader(filename=vcf_path)
    for i, vr in enumerate(reader):
        if limit is not None and i >= limit:
            break

        if not vr.ALT or vr.ALT[0] is None:
            continue

        if len(vr.ALT) > 1:
            raise ValueError("VCF contains multi-allelics: " + vr.line)

        record = {
            'chrom': vr.CHROM,
            'pos': vr.POS,
            'id': vr.ID,
            'ref': vr.REF,
            'alt': ",".join(map(str, vr.ALT)),
        }

        for key, value in vr.INFO.items():
            if key in ("OLD_MULTIALLELIC", "OLD_VARIANT", ):
                continue
            if isinstance(value, list):
                record[key] = value[0]
            else:
                record[key] = value

        yield record
