import argparse
import collections
import gzip
import os
import vcf

p = argparse.ArgumentParser()
p.add_argument("input_vcf_path")
p.add_argument("output_vcf_path")
args = p.parse_args()


input_filename_prefix = os.path.basename(args.input_vcf_path).replace(".vcf", "").replace(".gz", "")

print(f" Input: {args.input_vcf_path}")
print(f"Output: {args.output_vcf_path}")

with open(args.output_vcf_path, "w") as f:
    reader = vcf.VCFReader(filename=args.input_vcf_path)

    writer = vcf.VCFWriter(f, reader)

    counters = collections.defaultdict(int)
    for vr in reader:
        CallData = next(iter(reader._format_cache.values()))

        counters["total variants"] += 1
        if len(vr.ALT) > 1:
            p.error("VCF cannot have multi-allelics")

        ref = vr.REF
        if not ref or not vr.ALT or not vr.ALT[0]:
            continue

        alt = vr.ALT[0]
        samples = vr.samples

        repeat_unit = vr.INFO.get('RU', '').upper()

        id_field = ""
        if repeat_unit:
            id_field += "RU{}:{}:".format(len(repeat_unit), repeat_unit)

        if len(ref) != len(alt):
            if len(ref) > len(alt):
                id_field += "DEL:"
            elif len(ref) < len(alt):
                id_field += "INS:"

            id_field += "{}=>{}".format(len(ref)-1, len(alt)-1)

        record = vcf.model._Record(
            vr.CHROM,
            vr.POS,
            id_field,
            vr.REF,
            vr.ALT,
            vr.QUAL,
            vr.FILTER,
            vr.INFO,
            vr.FORMAT,
            sample_indexes={c.sample: c for c in samples},
            samples=samples)

        writer.write_record(record)


for key, value in sorted(counters.items()):
    print(f"{value:10d} {key}")

print("Done writing to " + args.output_vcf_path)
