import argparse
import collections
import os
import tqdm
import vcf

p = argparse.ArgumentParser()
p.add_argument("-s", "--replace-sample-name-with-filename", action="store_true")

g = p.add_mutually_exclusive_group(required=True)
g.add_argument("--gangstr", action="store_true")
g.add_argument("--expansion-hunter", action="store_true")

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
    for vr in tqdm.tqdm(reader, unit=" records"):
        CallData = next(iter(reader._format_cache.values()))

        counters["input lines"] += 1
        if not vr.ALT or vr.ALT[0] is None:
            continue

        counters["output lines"] += 1

        alt_alleles = vr.ALT
        samples = vr.samples

        if args.replace_sample_name_with_filename:
            updated_samples = []
            for sample in samples:
                updated_samples.append(
                    vcf.model._Call(sample.site, input_filename_prefix, sample.data))

            samples = updated_samples


        ref_allele = vr.REF
        if args.expansion_hunter:
            repeat_unit = vr.INFO['RU'] = vr.INFO['RU'].upper()
            ref_allele += int(vr.INFO['REF']) * repeat_unit

        # Both ExpansionHunter and GangSTR generate variants that have 2 identical alt alleles.
        # This collapses the identical alt-alleles into 1.
        if len(vr.ALT) == 2 and str(vr.ALT[0]) == str(vr.ALT[1]):
            alt_alleles = [vr.ALT[0]]
            counters["variants with alt1 == alt2"] += 1

            updated_samples = []
            for sample in samples:
                # replace 1/2 genotypes with 1/1, etc. (any  '2' genotype gets replace with a '1')
                new_sample_data = sample.data._asdict()
                new_sample_data['GT'] = new_sample_data['GT'].replace('2', '1')
                if args.replace_sample_name_with_filename:
                    sample_name = os.path.basename(args.input_vcf_path).replace(".vcf", "").replace(".gz", "")
                else:
                    sample_name = sample.sample

                updated_samples.append(
                    vcf.model._Call(sample.site, sample_name, CallData(*new_sample_data.values())))

            samples = updated_samples


        if args.expansion_hunter:
            repeat_unit = vr.INFO['RU'] = vr.INFO['RU'].upper()
            updated_alt_alleles = []
            for alt in alt_alleles:
                num_repeats_in_alt_allele = int(str(alt).strip('<STR>'))
                alt_allele = str(vr.REF) + repeat_unit * num_repeats_in_alt_allele
                if ref_allele == alt_allele:
                    counters["variants with ref_allele == alt_allele"] += 1

                updated_alt_alleles.append(vcf.model._Substitution(alt_allele))

            alt_alleles = updated_alt_alleles


        # rename END key to STR_END because END screws up how IGV display these variants.
        new_INFO = dict(vr.INFO)
        if 'END' in new_INFO:
            new_INFO['STR_END'] = new_INFO['END']
            del new_INFO['END']

        record = vcf.model._Record(
            vr.CHROM,
            vr.POS,
            vr.ID,
            ref_allele,
            alt_alleles,
            vr.QUAL,
            vr.FILTER,
            new_INFO,
            vr.FORMAT,
            sample_indexes={c.sample: c for c in samples},
            samples=samples)

        writer.write_record(record)


for key, value in sorted(counters.items()):
    print(f"{value:10d} {key}")

print("Done writing to " + args.output_vcf_path)