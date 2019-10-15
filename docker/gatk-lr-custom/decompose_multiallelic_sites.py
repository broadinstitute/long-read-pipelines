#!/usr/bin/env python2

# This script is copied, except these three lines, bit by bit from
# https://github.com/PacificBiosciences/hg002-ccs/blob/master/smallvariants/decompose_multiallelic_sites.py
# on 2019-09-27. We do not claim authorship in any means.

# UPDATE:
# It was found later that the line "annotation, values = infofield.split('=')"
# in function "decompose_INFO(INFO, position)" didn't precheck of when no "=" is present,
# so we added one before that.

"""Decompose multiallelic variant sites and statistics into biallelic sites.
"""

import os
import sys
import argparse


# FORMAT fields requiring recalculation
FORMAT_remove = ['PL']


def decompose_INFO(INFO, position):
    """Keep only data relevant to position.
    """
    newinfofields = []
    infofields = INFO.split(';')
    for infofield in infofields:
        # handling FLAG type INFO fields
        if "=" not in infofield:
            newinfofields.append(infofield)
            continue
        annotation, values = infofield.split('=')
        valuelist = values.split(',')
        if len(valuelist) == 2:
            # Split annotations that refer to one allele.
            # example: AC=1,1
            newinfofields.append('='.join([annotation, valuelist[position]]))
        elif annotation == 'AN':
            # For now, toss AN until it turns out we need it.
            pass
        elif len(valuelist) == 1:
            # Keep annotations that refer to the entire site.
            # example: QD=30.40
            newinfofields.append(infofield)
    return ';'.join(newinfofields)


def decompose_FORMAT(FORMAT):
    """Remove fields that require recalculation.
    Remove fields in the FORMAT_remove list.
    >>> decompose_FORMAT('GT:AD:DP:GQ:PL')
    'GT:AD:DP:GQ'
    """
    formatfields = FORMAT.split(':')
    return ':'.join([x for x in formatfields if x not in FORMAT_remove])


def decompose_GT(GT, FORMAT, position):
    """Transform multiallelic GT to biallelic GT
    Remove fields that require calculation, and take relevant subset from
    fields with data for multiple ALTs.
    >>> decompose_GT('1/2:1,8,12:21:99:766,319,340,191,0,201',
                     'GT:AD:DP:GQ:PL', 0)
    '0/1:1,8:21:99'
    """
    newgtfields = []
    gtfields = GT.split(':')
    formatfields = FORMAT.split(':')
    for field in zip(formatfields, gtfields):
        if field[0] in FORMAT_remove:
            pass
        elif field[0] == 'GT':
            newgtfields.append('0/1')
        elif field[0] == 'AD':
            values = field[1].split(',')
            newgtfields.append(','.join([values[0], values[position+1]]))
        else:
            newgtfields.append(field[1])
    return ':'.join(newgtfields)


def decompose_multiallelic(row):
    """Split, decompose, reassemble, and return rows.
    #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT GT
    Returns a list of new rows, one per ALT allele.
    """
    outrows = []
    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, GT = row.split()
    for position, alt in enumerate(ALT.split(',')):
         info = decompose_INFO(INFO, position)
         fmat = decompose_FORMAT(FORMAT)
         gt = decompose_GT(GT, FORMAT, position)
         outrows.append('\t'.join([CHROM, POS, ID,
                                   REF, alt, QUAL,
                                   FILTER, info, fmat, gt]))
    return outrows


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('invcf', help="Input file", type=argparse.FileType('r'))
    parser.add_argument('-o', '--outfile', help="Output file",
                    default=sys.stdout, type=argparse.FileType('w'))

    args = parser.parse_args(arguments)

    for row in args.invcf:
        # write all headers
        if row.startswith('#'):
            args.outfile.write(row)
        # decompose sites and write out new rows
        else:
            for outrow in decompose_multiallelic(row):
                args.outfile.write(outrow + '\n')


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
