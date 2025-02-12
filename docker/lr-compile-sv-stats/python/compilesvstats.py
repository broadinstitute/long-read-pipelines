import os
import subprocess
import argparse


def main():
    parser = argparse.ArgumentParser(description="""This script compiles the SV stats for a list of samples and a list of callers.""")
    parser.add_argument('--base_dir', required=True, help='The directory containing the sample files')
    parser.add_argument('--callers', nargs='+', required=True, help='The callers to compile stats for. (e.g.  --callers DELLY Strelka Manta)')
    parser.add_argument('--sample_file', required=True, help='The file containing the list of samples to compile stats for')
    parser.add_argument('--svtypes', nargs='+', default=["ALL", "DEL", "DUP", "CNV", "INS", "INV", "OTH"], help='The SV types to compile stats for. (e.g.  --svtypes ALL DEL DUP CNV INS INV OTH)')

    args = parser.parse_args()
    base_dir: str = args.base_dir
    callers: list = args.callers
    sample_file: str = args.sample_file
    svtypes: list = args.svtypes

    # resolve absolute paths
    resolved_base_dir = os.path.abspath(base_dir)
    resolved_sample_file = os.path.abspath(sample_file)

    # Check if paths exist
    if not os.path.isdir(resolved_base_dir):
        raise FileNotFoundError(f"base_dir: {resolved_base_dir} does not exist")
    if not os.path.isfile(resolved_sample_file):
        raise FileNotFoundError(f"sample_file: {resolved_sample_file} does not exist")


    basedir = resolved_base_dir

    with open(resolved_sample_file, 'r') as samplefile:
        samples = [line.strip() for line in samplefile]

        for caller in callers:
            with open("%s/%s_all_sample_stats" % (os.getcwd(), caller), 'w') as outfile:
                outfile.write("sample\t%s\n" % '\t'.join(svtypes))

                compile_stats(caller, svtypes, samples, basedir, outfile)


def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a)


def compile_stats(caller, svtypes, samples, basedir, outfile):
    for sample_name in samples:
        # Count lines in the file
        sample_svlen_file = os.path.join(basedir, f"{sample_name}.{caller}.txt")
        with open(sample_svlen_file, 'r') as file:
            ALL = sum(1 for _ in file)

        # Count SVs by type
        cut_sort_uniq_command = "cut -f1 {} | sort | uniq -c".format(sample_svlen_file)
        counts_by_sv = subprocess.check_output(cut_sort_uniq_command, shell=True)
        counts_by_sv_clean = counts_by_sv.decode().split()

        SVs = {}
        for svtype in svtypes:
            SVs[svtype] = 0
        SVs['ALL'] = ALL

        # Process counts_by_sv
        for i in range(0, len(counts_by_sv_clean), 2):
            num = int(counts_by_sv_clean[i])
            SV = counts_by_sv_clean[i + 1].upper()
            if SV in SVs:
                SVs[SV] = num
            else:
                SVs['OTH'] += num

        outfile.write("%s" % sample_name)

        for svtype in svtypes:
            outfile.write("\t%d" % SVs[svtype])
        outfile.write("\n")


if __name__ == "__main__":
    main()
