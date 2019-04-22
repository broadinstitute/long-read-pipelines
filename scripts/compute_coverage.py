#!/usr/bin/env python

# from Aaron Quinlin (https://github.com/arq5x/codachrom/blob/master/windowed_coverage.py)

import sys
import pybedtools as pbt
import argparse
import numpy as np
from collections import defaultdict


def make_windows(args):
    windows = pbt.BedTool().window_maker(
                                         w=args.window_size,
                                         s=args.step_size)

    tmp_file = "windows.temp.bed"
    x = windows.moveto(tmp_file)
    return tmp_file


def get_gc(args, windows_file):
    nuc_content = pbt.BedTool().nucleotide_content(fi=args.fasta_file, bed=windows_file)
    tmp_file = "nuc.temp.bed"
    x = nuc_content.moveto(tmp_file)
    return tmp_file


def get_rawcov(args, windows_file):
    windows = pbt.BedTool(windows_file)
    cov = pbt.BedTool(args.bam).coverage(windows, counts=True)
    return cov.sort()


class CoverageInterval(object):

    def __init__(self, ivl):
        self.ivl = ivl

    def set_gc(self, gc):
        self.gc = gc


def main():
    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(prog='chromcopy')
    parser.add_argument("-v", "--version", help="Installed chromcopy version",
                        action="version")
    parser.add_argument("-b", dest='bam', help="The BAM file for which to compute coverage.")
    parser.add_argument("-w", dest='window_size', help="The window size to be used when computing windows.")
    parser.add_argument("-s", dest='step_size', help="The step size to be used when computing windows.", default=0)
    parser.add_argument("-f", dest='fasta_file', help="The FASTA file of the genome you are using.")

    args = parser.parse_args()

    if args.bam is None or args.window_size is None or args.fasta_file is None:
        sys.exit('EXITING. You must specify -b, -f, and -w.\n')

    sys.stderr.write('Making windows...\n')
    windows_fn = make_windows(args)

    sys.stderr.write('Computing GC content in each window...\n')
    nuc_content_fn = get_gc(args, windows_fn)

    sys.stderr.write('Tabulating coverage in each window...\n')
    raw_cov = get_rawcov(args, windows_fn)

    # join the GC and counts
    gc_values = [float(feature[4]) for feature in pbt.BedTool(nuc_content_fn)]
    cov_intervals = []
    for idx, feature in enumerate(raw_cov):
        ivl = CoverageInterval(feature)
        ivl.set_gc(round(gc_values[idx], 2))
        cov_intervals.append(ivl)

    gc_bias_fh = open(args.bam + ".gcbias.txt", 'w')
    counts_per_gc = defaultdict(list)
    for ivl in cov_intervals:
        depth = int(ivl.ivl.name)
        counts_per_gc[ivl.gc].append(depth)
        gc_bias_fh.write(str(ivl.gc) + "\t" + str(depth) + '\n')

    mean_per_gc = defaultdict(float)
    stdv_per_gc = defaultdict(float)
    for gc in counts_per_gc:
        mean_per_gc[gc] = np.average(counts_per_gc[gc])
        stdv_per_gc[gc] = np.std(counts_per_gc[gc])

    depth_fh = open(args.bam + ".depth.txt", 'w')
    for ivl in cov_intervals:
        gc = ivl.gc
        mean = mean_per_gc[gc]
        stdv = stdv_per_gc[gc]
        Z = float(int(ivl.ivl.name) - mean) / float(stdv + 1.0)
        depth_fh.write('\t'.join(
            [str(s) for s in [ivl.ivl.chrom, ivl.ivl.start, ivl.ivl.end, ivl.ivl.name, ivl.gc, mean, stdv, Z]]) + '\n')

    # cleaning up.
    sys.stderr.write('Cleaning up...\n')
    os.remove(windows_fn)
    os.remove(nuc_content_fn)


if __name__ == "__main__":
    main()
