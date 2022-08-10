import argparse
from typing import List
import pandas as pd
import numpy as np


def bin_and_freq(dyst_list:str) -> pd.Series:
    try:
        arr = dyst_list.split('-')
        beg = arr[0].strip()
        brr = arr[1].split('[')
        end = brr[0].strip()
        freq = brr[1].replace(']', '').strip()
        return pd.Series([float(beg), float(end), int(freq)], index = ['beg', 'end', 'freq'])
    except:
        print(dyst_list)
        raise


def construct_dataframe(prepped_dyst_output:str) -> pd.DataFrame:
    """
    Converts a prepped dyst output format into a dataframe with three columns,
    [left boundary of a bin, right boundary of a bin, frequency of the bin]
    :param prepped_dyst_output:
    :return:
    """
    with open(prepped_dyst_output, 'r') as inf:
        lines = [l.strip('\n') for l in inf.readlines()]
    df = pd.DataFrame([bin_and_freq(l) for l in lines if l])
    df['freq'] = df['freq'].astype(int)
    return df


def find_blocks_in_sequence(indices: List[int]) -> List[List[int]]:
    """
    Given a list of indices into a dataframe, produce blocks of contiguous indices.
    :param indices:
    :return:
    """
    assert all([i>=0 for i in indices]), "input value must be all non-negative"
    blocks = list()
    indices_shift_1 = indices[1:]
    indices_shift_1.append(indices[-1]+1)
    i = 0
    d = 1
    for a, b in zip(indices, indices_shift_1):
        if b-a!=1:
            blocks.append(indices[i:(i+d)])
            i+=d
            d=0
        d+=1
    blocks.append(indices[i:(i+d)])
    return blocks


def heuristic_trough_to_peak_ratio(formatted_df:pd.DataFrame, ratio:float) -> List[int]:
    """
    The heuristic part of this algorithm.

    It works by finding the bins whose frequency is higher than _ratio_ * the globally highest frequency.
    This ratio is a heuristic number.
    :param formatted_df:
    :param ratio:
    :return:
    """
    sequence = list(formatted_df.loc[formatted_df['freq'] > round(max(formatted_df['freq'])*ratio), :].index)
    return sequence


def find_peaks_indices(formatted_df:pd.DataFrame, blocks_of_contiguous_indices:List[List[int]]) -> list:
    """
    Given the blocks of indices, where each block is contiguous within themselves,
    locate the index into the dataframe where the peak sits.
    :param formatted_df:
    :param blocks_of_contiguous_indices:
    :return:
    """
    indices_of_peaks = list()
    for bl in blocks_of_contiguous_indices:
        tdf = formatted_df.iloc[bl, :]
        idx = tdf['freq'].idxmax()
        indices_of_peaks.append(idx)
    return indices_of_peaks


def bootstrap_peak_indices(formatted_df:pd.DataFrame, heuristic_ratio_lower_bound:float, heuristic_ratio_upper_bound:float):

    assert 0 < heuristic_ratio_lower_bound < 1, \
        f"lower bound ({heuristic_ratio_lower_bound}) must be between 0 and 1."
    assert 0 < heuristic_ratio_upper_bound < 1, \
        f"lower bound ({heuristic_ratio_upper_bound}) must be between 0 and 1."
    assert heuristic_ratio_lower_bound < heuristic_ratio_upper_bound, \
        f"lower bound ({heuristic_ratio_lower_bound}) must be lower than upper bound ({heuristic_ratio_upper_bound})."

    result = list()
    # bootstrap
    for f in np.arange(heuristic_ratio_lower_bound, heuristic_ratio_upper_bound, 0.1):
        sequence = heuristic_trough_to_peak_ratio(formatted_df, f)
        blocks = find_blocks_in_sequence(sequence)
        peaks = find_peaks_indices(formatted_df, blocks)
        print(peaks)
        result.extend(peaks)
    return sorted(list(set(result)))


def qc_check_peak_indices(formatted_df:pd.DataFrame, bootstrapped_peaks: list):
    assert 0 < len(bootstrapped_peaks), "Please call me only when you've found peaks"

    if 1 == len(bootstrapped_peaks):
        return bootstrapped_peaks

    bootstrapped_peaks_shifted_right = bootstrapped_peaks[1:]
    bootstrapped_peaks_drop_last = bootstrapped_peaks[:-1]

    qc_pass_list = list()
    for a, b in zip(bootstrapped_peaks_drop_last, bootstrapped_peaks_shifted_right):
        x = formatted_df.iloc[a, :]['freq']
        y = formatted_df.iloc[b, :]['freq']
        m = min(x, y)
        if any([f < m for f in formatted_df.iloc[a+1:b, :]['freq']]):
            qc_pass_list.extend([a, b])
    return sorted(list(set(qc_pass_list)))


def locate_peak_bin_value(formatted_df:pd.DataFrame, peak_indices:List[int]) -> list:
    peaks = list()
    for idx in peak_indices:
        peaks.append(round((formatted_df.iloc[idx, :]['beg'] + formatted_df.iloc[idx, :]['end'])/2))
    return peaks


######################################################################
def main():
    parser = argparse.ArgumentParser(description='Find the first few peaks in a (prepped) dyst histogram output',
                                     prog='find_peaks')
    parser.add_argument('-i', '--input', type=str, help="prepped dyst histogram output (without the bars)")
    parser.add_argument('-o', '--output', type=str, help="path to output peak values (flat file)")
    args = parser.parse_args()

    prepped_dyst_output = args.input
    peaks_output = args.output
    df = construct_dataframe(prepped_dyst_output)
    peak_indices = bootstrap_peak_indices(df, 0.2, 0.8)
    pass_qc_peak_indices = qc_check_peak_indices(df, peak_indices)
    peak_values = locate_peak_bin_value(df, pass_qc_peak_indices)
    with open(peaks_output, 'w') as outf:
        [outf.write(f"{p}\n") for p in peak_values]


######################################################################
if __name__ == "__main__":
    main()

