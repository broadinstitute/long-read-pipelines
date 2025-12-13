import argparse
import gzip
import json
import logging

import numpy as np
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process coverage values.")
    parser.add_argument("mosdepth_regions", type=str,
                        help="Mosdepth output file containing coverage values.")
    parser.add_argument("--cov_col", type=int, default=4,
                        help="Column holding the coverage values.")
    parser.add_argument("--output_prefix", type=str,
                        help="Prefix to use for output files.")
    parser.add_argument("--round", type=int, default=2, help="Rounding precision.")
    parser.add_argument("--debug", action="store_true", help="Debug mode.")
    return parser.parse_args()


def calculate_summary_statistics(
        df: pd.DataFrame, cov_col: int, round_precision: int
) -> dict:
    """
    Calculate summary statistics for the coverage values.
    @param df: Dataframe containing coverage values
    @param cov_col: Column containing coverage values
    @param round_precision: Rounding precision
    @return:
    """

    cov_data: pd.Series = df.iloc[:, cov_col - 1]
    mean_val: float = cov_data.mean()
    mad_val: float = np.mean(np.abs(cov_data - mean_val))

    statistics = {
        "cov_mean": round(mean_val, round_precision),
        "cov_q1": round(cov_data.quantile(0.25), round_precision),
        "cov_median": round(cov_data.median(), round_precision),
        "cov_q3": round(cov_data.quantile(0.75), round_precision),
        "cov_iqr": round(cov_data.quantile(0.75) - cov_data.quantile(0.25), round_precision),
        "cov_stdev": round(cov_data.std(), round_precision),
        "cov_mad": round(mad_val, round_precision),
        "cov_percent_above_4x": percentage_greater_than_4x(df, cov_col, round_precision),
        "cov_evenness_score": calculate_evenness_score(df, cov_col, round_precision)
    }

    # Replace Nan values with null
    for key, value in statistics.items():
        if pd.isna(value):
            statistics[key] = "null"

    return statistics


def percentage_greater_than_4x(
        df: pd.DataFrame, cov_col: int, round_precision: int
) -> float:
    """
    Calculate the percentage of coverage values greater than 4x.
    @param df: Dataframe containing coverage values
    @param cov_col: Column containing coverage values
    @param round_precision: Rounding precision
    @return:
    """
    logging.debug("Calculating percentage of coverage values greater than 4x")
    total_bases = len(df)
    bases_above_4x = df[df.iloc[:, cov_col - 1] > 4].shape[0]
    percent_above_4x = round(bases_above_4x / total_bases, round_precision)

    logging.info("Percentage of coverage values greater than 4x: %s", percent_above_4x)

    return percent_above_4x


def calculate_evenness_score(
        df: pd.DataFrame, cov_col: int, round_precision: int) -> float:
    """
    Calculate the evenness score.
    Konrad Oexle, Journal of Human Genetics 2016, Evaluation of the evenness score in NGS.
    https://www.nature.com/articles/jhg201621

    @param df: Dataframe containing coverage values
    @param cov_col: Column containing coverage values
    @param round_precision: Rounding precision
    @return:
    """
    logging.debug("Calculating evenness score")
    mean_coverage = df.iloc[:, cov_col - 1].mean()
    # Get the coverage values that are less than or equal to the mean coverage
    d2 = df[df.iloc[:, cov_col - 1] <= mean_coverage].iloc[:, cov_col - 1].tolist()
    # count of coverage values that are less than or equal to the mean coverage
    d2_count = len(d2)
    # sum of coverage values that are less than or equal to the mean coverage
    d2_sum = sum(d2)
    # total number of coverages
    coverage_count = len(df)

    logging.debug("Mean coverage: %s", mean_coverage)
    logging.debug("D2 count: %s", d2_count)
    logging.debug("D2 sum: %s", d2_sum)
    logging.debug("Coverage count: %s", coverage_count)

    if mean_coverage != 0:
        evenness_score = round(1 - (d2_count - (d2_sum / mean_coverage)) / coverage_count, round_precision)
    else:
        logging.warning("Mean coverage is zero. Evenness score will be set to null.")
        evenness_score = "null"

    logging.info("Evenness score: %s", evenness_score)

    return evenness_score


def open_file(file_path: str) -> pd.DataFrame:
    """
    Open a file and return a dataframe
    @param file_path: Path to the file
    @return:
    """

    if file_path.endswith(".gz"):
        with gzip.open(file_path, "rt") as f:
            df = pd.read_csv(f, sep="\t", header=None)
    else:
        with open(file_path, "rt") as f:
            df = pd.read_csv(f, sep="\t", header=None)

    return df


def main():
    args = parse_arguments()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.info("Arguments: %s", args)
    logging.info("Calculating coverage statistics")

    prefix = args.output_prefix if args.output_prefix else args.mosdepth_regions.split(".")[0]

    df: pd.DataFrame = open_file(args.mosdepth_regions)
    logging.info("Opened file: %s", args.mosdepth_regions)

    logging.debug("Dataframe shape: %s", df.shape)
    logging.debug("Dataframe columns: %s", df.columns)
    logging.debug(f"Dataframe Head\n{df.head()}")

    statistics = calculate_summary_statistics(df, args.cov_col, args.round)

    logging.info("Summary statistics: %s", statistics)

    summary_file = f"{prefix}.cov_stat_summary.json"

    logging.info("Writing summary statistics to file: %s", summary_file)
    with open(summary_file, "w") as f:
        json.dump(statistics, f)
        f.write("\n")


if __name__ == "__main__":
    main()
