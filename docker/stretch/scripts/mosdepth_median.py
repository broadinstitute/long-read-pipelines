#!/usr/bin/env python
"""Estimate the median coverage from the mosdepth.global.dist.txt cumulative coverage output
"""

import argparse
import sys
import pandas as pd
import numpy as np

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description='Estimate allele lengths and find outliers at STR loci.')
    parser.add_argument(
        'INFILE', type=str, 
        help='mosdepth.global.dist.txt cumulative coverage output file')
    parser.add_argument(
        '--out', type=str, required=False,
        help='Output file name. Defaults to stdout.')
    return parser.parse_args()

def parse_dist(filename):
    """Parse mosdepth.global.dist.txt to pandas df"""
    try:
        data = pd.read_table(filename, delim_whitespace = True, 
            names = ['chr', 'coverage', 'prop_bases_covered'])
    except pd.io.common.EmptyDataError:
        sys.exit('ERROR: file {0} was empty.\n'.format(filename))
    return(data)

def find_nearest(array, value):
    """https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def main():
    args = parse_args()
    infile = args.INFILE
    outfile = args.out

    if outfile:
        outstream = open(outfile, 'w')
    else:
        outstream = sys.stdout

    dist_data = parse_dist(infile)
    # extract values from all chromosomes
    dist_data_total = dist_data.loc[dist_data['chr'] == 'total']
    # find the median coverage (the coverage where prop_bases_covered is closest to 0.5)
    closest_prop = find_nearest(dist_data_total['prop_bases_covered'], 0.5)
    median_cov = dist_data_total.loc[dist_data_total['prop_bases_covered'] == closest_prop]['coverage']

    outstream.write(str(median_cov.iat[0])+'\n')

if __name__ == '__main__':
    main()
