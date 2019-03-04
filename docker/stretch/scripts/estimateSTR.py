#!/usr/bin/env python
"""Estimate allele lengths and find outliers at STR loci
"""

import warnings
from pprint import pprint

with warnings.catch_warnings():
    warnings.simplefilter("ignore", FutureWarning)

    import argparse
    import sys
    import glob
    import os
    import re
    import numpy as np
    import statsmodels.api as sm
    from scipy.stats import norm
    from statsmodels.sandbox.stats.multicomp import multipletests
    from sklearn import linear_model
    import pandas as pd

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description='Estimate allele lengths and find outliers at STR loci.')
    parser.add_argument(
        '--locus_counts', type=str, nargs='+', required = True,
        help='.locus_counts files for all samples. Contains the number of reads assigned to each STR locus.')
    parser.add_argument(
        '--STR_counts', type=str, nargs='+', required = True,
        help='.STR_counts files for all samples. Contains the number of reads mapped to each STR decoy chromosome.')
    parser.add_argument(
        '--median_cov', type=str, nargs='+', required = True,
        help='.median_cov files for all samples. Text files containing median coverage.')
    parser.add_argument(
        '--out', type=str, default = '',
        help='Prefix for all output files (suffix will be STRs.tsv) (default: %(default)s)')
    parser.add_argument(
        '--model', type=str, default='STRcov.model.csv',
        help='Data to produce linear regression model (provided with STRetch) (default: %(default)s)')
    parser.add_argument(
        '--control', type=str, default='',
        help='Input file for median and standard deviation estimates at each locus from a set of control samples. This file can be produced by this script using the emit option. If this option is not set, all samples in the current batch will be used as controls by default.')
    parser.add_argument(
        '--emit', type=str, default='',
        help='Output file for median and standard deviation estimates at each locus (tsv).')
    return parser.parse_args()

def get_sample(fullpath):
    """Get the sample ID from the filename"""
    basename = os.path.basename(fullpath)
    return(basename.split('.')[0])

def parse_STRcov(filename):
    """Parse all STR coverage"""
    sample_id = get_sample(filename)
    try:
        cov_data = pd.read_table(filename, delim_whitespace = True,
                            names = ['chrom', 'start', 'end', 'decoycov'])
    except pd.io.common.EmptyDataError:
        sys.exit('ERROR: file {0} was empty.\n'.format(filename))
    cov_data['sample'] = sample_id
    cov_data['repeatunit'] = [x.split('-')[1] for x in cov_data['chrom']]
    cov_data = cov_data[['sample', 'repeatunit', 'decoycov']]
    return(cov_data)

def parse_locuscov(filename):
    """Parse locuscoverage data produced by identify_locus.py"""
    sample_id = get_sample(filename)
    try:
        locuscov_data = pd.read_table(filename, delim_whitespace = True)
    except pd.io.common.EmptyDataError:
        sys.exit('ERROR: file {0} was empty.\n'.format(filename))
    if locuscov_data.shape[0] == 0: # Check for file with only header
        sys.exit('ERROR: file {0} contained 0 loci.\n'.format(filename))

    locuscov_data['sample'] = sample_id
    locuscov_data['locus'] = ['{0}-{1}-{2}'.format(locuscov_data['STR_chr'][i],
                            locuscov_data['STR_start'][i], locuscov_data['STR_stop'][i]) for i in range(len(locuscov_data.index-1))]
    locuscov_data['repeatunit'] = locuscov_data['motif']
    locuscov_data['locuscoverage'] = locuscov_data['count']
    locuscov_data = locuscov_data[['sample', 'locus', 'repeatunit', 'reflen', 'locuscoverage']]
    return(locuscov_data)

def parse_genomecov(filename):
    """Parse median genome coverage from covmed output.
        Assumes median coverage is the top left value in the text file."""
    sample_id = get_sample(filename)
    try:
        mediancov = pd.read_table(filename, delim_whitespace = True, header = None).iloc[0,0]
    except pd.io.common.EmptyDataError:
        sys.exit('ERROR: file {0} was empty.\n'.format(filename))
    if mediancov < 1:
        sys.exit('ERROR: Median coverage in file {0} was {1}.\nSuch a low value could indicate median coverage was not correctly calculated,\nfor example an incorrect target region was specified or the WGS pipeline was used for exome data.'.format(filename, mediancov))
    genomecov_data = pd.DataFrame({'sample': [sample_id], 'genomecov': [mediancov]})
    return(genomecov_data)

def parse_controls(control_file):
    """Parse control file with columns locus, median and standard deviation"""

    control_estimates = pd.read_table(control_file, index_col=0)

    # Allow for old style column headings, but change to mu and sd.
    if control_estimates.columns[0] in ['mu', 'median'] and control_estimates.columns[1] in ['sd', 'SD']:
        colnames = list(control_estimates.columns)
        colnames[0:2] = ['mu', 'sd']
        control_estimates.columns = colnames
    else:
        raise ValueError(''.join(["The column names in the control file ",
        "don't look right, expecting columns named median, SD ",
        "or mu, sd. Column names are ", str(list(control_estimates.columns)),
        ". Check the file: ", control_file]))
    return(control_estimates)

#from statsmodels import robust
# If using mad below

def hubers_est(x):
    """Emit Huber's M-estimator median and SD estimates.
    If Huber's fails, emit standard median and NA for sd"""
    huber50 = sm.robust.scale.Huber(maxiter=50)

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)

        try:
            mu, s = huber50(np.array(x))
        except (ValueError, RuntimeWarning):
            mu = np.median(x)
            s = np.nan
            #s = robust.mad(x)
            #XXX working on this - replace s with mad when hubers est fails?

    return pd.Series({'mu': mu, 'sd': np.sqrt(s)})

def z_score(x, df):
    """Calculate a z score for each x value, using estimates from a pandas data
    frame with the columns 'mu' and 'sd' and index coressponding to the x values"""
    z = (x.transpose() - df['mu'])/df['sd']
    return z.transpose()

def p_adj_bh(x):
    '''Adjust p values using Benjamini/Hochberg method'''
    return multipletests(x, method='fdr_bh', returnsorted = False)[1]

def main():
    # Parse command line arguments
    args = parse_args()

    print("===> estimate_STR.py input args:")
    pprint(args.__dict__)

    base_filename = args.out
    STRcov_model_csv = args.model
    emit_file = args.emit
    control_file = args.control
    locuscov_files = args.locus_counts
    STRcov_files = args.STR_counts
    genomecov_files = args.median_cov
    results_suffix = 'STRs.tsv'

    # Check files exist for all samples
    locuscov_ids = set([get_sample(f) for f in locuscov_files])
    STRcov_ids = set([get_sample(f) for f in STRcov_files])
    genomecov_ids = set([get_sample(f) for f in genomecov_files])
    if not (locuscov_ids == STRcov_ids == genomecov_ids):
        all_samples = locuscov_ids | STRcov_ids | genomecov_ids
        missing_samples = (all_samples - locuscov_ids) | (all_samples - STRcov_ids) | (all_samples - genomecov_ids)
        sys.exit("ERROR: One or more files are missing for sample(s): " + ' '.join(missing_samples))
    
    sys.stderr.write('Processing {0} samples\n'.format(len(locuscov_files)))

    if len(locuscov_files) < 2 and control_file == '':
        sys.stderr.write('WARNING: Only 1 sample and no control file provided, so outlier scores and p-values will not be generated.\n')

    # Parse input data
    locuscov_data = pd.concat( (parse_locuscov(f) for f in locuscov_files), ignore_index = True)
    STRcov_data = pd.concat( (parse_STRcov(f) for f in STRcov_files), ignore_index = True)
    genomecov_data = pd.concat( (parse_genomecov(f) for f in genomecov_files), ignore_index = True)

    # Check for multiple rows with the same sample/locus combination
    crosstable = pd.crosstab(locuscov_data['locus'], locuscov_data['sample'])
    ismultiloci = crosstable.apply(lambda row: any(row > 1), axis=1)
    multiloci = ismultiloci[ismultiloci == True].index.values
    if len(multiloci) > 0:
        sys.exit('''
    The locus count input data contains multiple rows with the same sample/locus combination.
    This is usually caused by two loci at the same position in the STR annotation bed file.
    Check these loci:
    ''' + ' '.join(multiloci))
    
#    # Check for different reflen for the same locus
#    grouped = locuscov_data.groupby('locus')
#    reflenloci = []
#    for locus, group in grouped:
#        if len(set(group['reflen'])) > 1:
#            #reflenloci.append(name)
#            # If different, replace with the smallest
#            locuscov_data.loc[locuscov_data['locus'] == locus,'reflen'] = np.repeat(min(group['reflen']), len(group['reflen']))
#    if len(reflenloci) > 0:
#         sys.exit('''
#    The locus count input data contains the same locus with different reflens.
#    This may be caused by an error in the STR annotation bed file.
#    Check these loci:
#    ''' + ' '.join(reflenloci)) + '''
#    The locus count input data contains the same locus with different reflens.
#    This may be caused by an error in the STR annotation bed file.
#    Check the above loci'''

    #locuscov_data['reflen'] = np.repeat(1, len(locuscov_data['reflen']))

    # Fill zeros in locuscov
    locuscov_wide = locuscov_data.pivot(index='locus', columns='sample', values='locuscoverage').fillna(0)
    locuscov_wide['locus'] = locuscov_wide.index

    sample_cols = list(set(locuscov_data['sample']))
    locuscov_long = pd.melt(locuscov_wide, id_vars = 'locus',
                            value_vars = sample_cols, value_name = 'locuscoverage',
                            var_name = 'sample')

    # Add locus info back in
    locuscov_data = pd.merge(locuscov_long, locuscov_data[['locus', 'repeatunit', 'reflen']].drop_duplicates(), how='left')

    # Normalise STR coverage by median coverage
    factor = 100

    STRcov_data = pd.merge(STRcov_data, genomecov_data)
    #STRcov_data['decoycov_norm'] = factor * (STRcov_data['decoycov'] + 1) / STRcov_data['genomecov']
    #STRcov_data['decoycov_log'] = np.log2(STRcov_data['decoycov_norm'])
    #XXX combines the previous two lines into one. Keeping commented out in case coverage_norm is required later
    STRcov_data['decoycov_log'] = np.log2(factor * (STRcov_data['decoycov'] + 1) / STRcov_data['genomecov'])

    print("==> STRcov_data shape: " + str(STRcov_data.shape))

    # #XXX ***Not fully implemented***
    # # Look for repeat units where the number of reads mapping to the decoy can't be
    # # explained by those mapping to all loci with that repeat unit
    #
    # # Sum the counts over all loci for each repeat unit in each sample
    # locus_totals = locuscov_data.groupby(['sample', 'repeatunit'])['locuscoverage'].aggregate(np.sum)
    # locus_totals = pd.DataFrame(locus_totals).reset_index() # convert to DataFrame and make indices into columns
    # # Calculate the difference between reads assigned to a decoy and the sum of
    # # all reads assigned to loci with that repeat unit
    # all_differences = pd.merge(STRcov_data, locus_totals, how='left')
    # all_differences['difference'] = all_differences['decoycov'] - all_differences['locuscoverage']
    # # Normalise differences by median coverage and take the log2
    # all_differences['difference_log'] = np.log2(factor * (all_differences['difference'] + 1) / all_differences['genomecov'])
    #
    # locus_totals = pd.merge(locus_totals, STRcov_data)
    #
    # # Assign decoy counts to each locus, based on what proportion of the counts for that repeat unit they already have

    locus_totals = pd.merge(locuscov_data, STRcov_data, how = 'left')

    locus_totals['total_assigned'] = locus_totals['locuscoverage'] #XXX remove this line if implementing the above
    locus_totals['total_assigned_log'] = np.log2(factor * (locus_totals['total_assigned'] + 1) / locus_totals['genomecov'])
    
    # For each locus, calculate if that sample is an outlier relative to others
    total_assigned_wide = locus_totals.pivot(index='locus', columns='sample', values='total_assigned_log')

    # Calculate values for if there were zero reads at a locus in all samples
    null_locus_counts = np.log2(factor * (0 + 1) / genomecov_data['genomecov'])
    sample_names = genomecov_data['sample']
    null_locus_counts.index = sample_names
    # Add a null locus that has 0 reads for all individuals
    # (so just uses coverage)
    null_locus_counts_est = hubers_est(null_locus_counts)

    # Calculate a z scores using median and SD estimates from the current set
    # of samples

    # Use Huber's M-estimator to calculate median and SD across all samples
    # for each locus
    sample_estimates = total_assigned_wide.apply(hubers_est, axis=1)
    # Where sd is NA, replace with the minimum non-zero sd from all loci
    min_sd = np.min(sample_estimates['sd'][sample_estimates['sd'] > 0])
    sample_estimates['sd'].fillna(min_sd, inplace=True)

    # if sd is 0, replace with min_sd #XXX is this sensible?
    #if null_locus_counts_est['sd'] == 0 or np.isnan(null_locus_counts_est['sd']):
    if null_locus_counts_est['sd'] == 0:
        null_locus_counts_est['sd'] = min_sd

    # Save median and SD of all loci to file if requested (for use as a
    # control set for future data sets)
    if emit_file != '':

        sample_estimates.loc['null_locus_counts'] = null_locus_counts_est

        n = len(total_assigned_wide.columns)
        sample_estimates['n'] = n

        sample_estimates.to_csv(emit_file, sep= '\t')

    # Calculate z scores using median and SD estimates per locus from a
    # provided control set
    if control_file != '':
        # Parse control file
        control_estimates = parse_controls(control_file)
        # Get a list of all loci in the control file but not the sample data
        control_loci_df = control_estimates.iloc[control_estimates.index != 'null_locus_counts']
        control_loci = [x for x in control_loci_df.index if x not in total_assigned_wide.index]

        # Extract and order just those control estimates appearing in the current data
        mu_sd_estimates = control_estimates.reindex(total_assigned_wide.index)
        # Fill NaNs with null_locus_counts values
        mu_sd_estimates.fillna(control_estimates.loc['null_locus_counts'],
                                inplace=True)
    else:
        # Extract and order estimates to match the current data
        mu_sd_estimates = sample_estimates.reindex(total_assigned_wide.index)

    # calculate z scores
    z = z_score(total_assigned_wide, mu_sd_estimates)

    # If a control file is given, effectively add zeros counts at all loci in 
    # controls but not in the samples. 
    # These extra rows will dissapear due to a later merge
    if control_file != '': 
        # Create a total_assigned_wide as if all loci have zero counts
        null_total_assigned_wide = pd.DataFrame(columns = sample_names, index = control_loci)
        null_total_assigned_wide.fillna(null_locus_counts, inplace = True)
        # Caculate z scores
        null_z = z_score(null_total_assigned_wide, 
                            control_estimates.reindex(null_total_assigned_wide.index))
        loci_with_counts = z.index
        z = z.append(null_z)

    if z.shape[0] == 1:
        ids = z.columns # save index order as data gets sorted
        # Calculate p values based on z scores (one sided)
        z_list = list(z.iloc[0])
        pvals = norm.sf(z_list) # no need to adjust p values if one locus

        # Merge pvals and z scores back into locus_totals
        p_z_df = pd.DataFrame({'sample': ids, 'p_adj': pvals, 'outlier': z_list})
        locus_totals = pd.merge(locus_totals, p_z_df)

    elif z.shape[0] > 1:
        # Calculate p values based on z scores (one sided)
        pvals = z.apply(lambda z_row: [norm.sf(x) for x in z_row], axis=1, result_type='broadcast') # apply to each row

        if pvals.isnull().values.all(): # Don't bother adjusting p values if all are null
            adj_pvals = pvals
        else:
            # Adjust p values using Benjamini/Hochberg method
            adj_pvals = pvals.apply(p_adj_bh, axis=0) # apply to each column
        
        # Merge pvals and z scores back into locus_totals
        adj_pvals['locus'] = adj_pvals.index
        pvals_long = pd.melt(adj_pvals, id_vars = 'locus',
                                value_vars = sample_cols, value_name = 'p_adj', var_name = 'sample')
        locus_totals = pd.merge(locus_totals, pvals_long)
        
        z['locus'] = z.index #important to do this only after p values calculated
        z_long = pd.melt(z, id_vars = 'locus',
                        value_vars = sample_cols, value_name = 'outlier', var_name = 'sample')
        locus_totals = pd.merge(locus_totals, z_long)

    elif z.shape[0] == 0:
        pass #XXX raise error. No input data!

    # Predict size (in bp) using the ATXN8 linear model (produced from data in
    # decoySTR_cov_sim_ATXN8_AGC.R)
    # Read in the raw data for this model from a file
    # Note: coverage_norm = (STR coverage/median coverage) * 100
    # allele2 is the length of the longer allele in bp inserted relative to ref
    STRcov_model = pd.read_csv(STRcov_model_csv)

    # Model is built from log2 data then converted back
    # (to reduce heteroscedasticity)

    # Create linear regression object
    regr = linear_model.LinearRegression()
    # Train the model
    # Reshape using X.reshape(-1, 1) if data has a single feature
    # or X.reshape(1, -1) if it contains a single sample.
    X_train = np.log2(STRcov_model['coverage_norm']).values.reshape(-1, 1)
    Y_train = np.log2(STRcov_model['allele2'])
    regr.fit(X_train, Y_train)
    # Make a prediction
    Y_pred = regr.predict(locus_totals['total_assigned_log'].values.reshape(-1, 1))
    predict = np.power(2, Y_pred)
    locus_totals['bpInsertion'] = predict

    # Get the estimated size in terms of repeat units (total, not relative to ref)
    repeatunit_lens = [len(x) for x in locus_totals['repeatunit']]
    locus_totals['repeatUnits'] = (locus_totals['bpInsertion']/repeatunit_lens) + locus_totals['reflen']

    # Split locus into 3 columns: chrom start end
    locuscols = pd.DataFrame([x.split('-') for x in locus_totals['locus']],
                        columns = ['chrom', 'start', 'end'])
    locus_totals = locus_totals.join(locuscols)

    # Specify output data columns
    write_data = locus_totals[['chrom', 'start', 'end',
                                    'sample', 'repeatunit', 'reflen',
                                    'locuscoverage',
                                    'outlier', 'p_adj',
                                    'bpInsertion', 'repeatUnits'
                                    ]]

    #sort by outlier score then estimated size (bpInsertion), both descending
    write_data = write_data.sort_values(['outlier', 'bpInsertion'], ascending=[False, False])
    #XXX check for duplicate rows?

    # Write individual files for each sample, remove rows where locuscoverage == 0
    samples = set(write_data['sample'])
    for sample in samples:
        sample_filename = base_filename + sample + '.' + results_suffix
        sample_df = write_data.loc[write_data['sample'] == sample]
        sample_df = sample_df.loc[sample_df['locuscoverage'] != 0.0]
        sample_df.to_csv(sample_filename, sep= '\t', index = False, na_rep='NaN')

    # Write all samples to a single file
    all_filename = base_filename + results_suffix
    write_data.to_csv(all_filename, sep= '\t', index = False, na_rep='NaN')

if __name__ == '__main__':
    main()
