#!/usr/bin/env python
# coding: utf-8

# Processes barcode data for P. falciparum samples
# Based on work from the following publication: https://doi.org/10.1093/pnasnexus/pgac187
# Orignal Author: Wesley Wong
# Modified by: Jonn Smith

import random
import copy
import itertools
import sys
import json
import argparse

import scipy
import pandas as pd
import numpy as np
import seaborn as sns
import numpy.polynomial.polynomial as poly
import statsmodels.api as sm
import networkx as nx
import matplotlib.pyplot as plt

from pandas import DataFrame, read_csv, read_excel
from collections import defaultdict, Counter
from matplotlib.patches import Patch
from patsy import dmatrices
from sklearn import linear_model
from numpyencoder import NumpyEncoder
from sklearn.neighbors import KernelDensity

def remove_duplicate_df(df, field="Haplotype"):
    duplicate_rows = pd.DataFrame.duplicated(df, field, keep="first")
    return df[~duplicate_rows]


def inverse_var_weight(p_array, var_array):
    var_weighted_num = np.nansum(np.divide(p_array, var_array, where=var_array != 0))
    var_weighted_denom = np.nansum(
        np.divide([1 for _ in var_array], var_array, where=var_array != 0)
    )
    weighted_mean = var_weighted_num / var_weighted_denom

    weighted_var = 1 / np.sum(1 / var_array)
    weighted_std = np.sqrt(weighted_var)
    weighted_ci = (
        weighted_mean - 1.96 * weighted_std,
        weighted_mean + 1.96 * weighted_std,
    )

    return weighted_mean, weighted_var, weighted_ci


def is_unique(s):
    a = s.to_numpy()  # s.values (pandas<0.24)
    return (a[0] == a).all()  # True = homozygous, #False = Heterozygous


def return_stats(s, convert_numpy=True):
    if convert_numpy:
        a = s.to_numpy()  # s.values (pandas<0.24)
    else:
        s = a
    missing_screen = a != "X"
    poly_screen = a != "N"

    if missing_screen.all() == False:
        return np.NaN  # missing data present in this comparison, skip
    elif poly_screen.all() == True:
        return 1.0 - (a[0] == a).all()  # 0 = homozygous, #1 = Heterozygous

    else:  # N present in one of the comparisons, skip for now
        return np.NaN


def quantify_het(barcode_comparison, f_true=0, axis=0):
    null_het = np.nansum(barcode_comparison, axis=axis) * (1 - f_true)
    return null_het


def quantify_total(barcode_comparison, axis=0):
    return np.sum(~np.isnan(barcode_comparison), axis=axis)


def run_fn_tests():
    test_het_calculator()
    test_unique()


def test_het_calculator():
    expected_het, expected_total = (12, 16)
    test_barcode1 = (
        6 * ["A"] + 6 * ["T"] + 6 * ["C"] + 6 * ["G"] + 6 * ["X"] + 6 * ["N"]
    )
    test_barcode2 = 6 * ["A", "T", "C", "G", "X", "N"]
    test_df = DataFrame([test_barcode1, test_barcode2])
    b = test_df.apply(return_stats, axis=0).to_list()
    heterozygotes = np.sum(quantify_het(b), axis=0)
    total = np.sum(quantify_total(b), axis=0)
    assert (
        heterozygotes == expected_het
    ), "Heterozygote sites should be {expected_het}, identified {x}".format(
        expected_het=expected_het, x=heterozygotes
    )
    assert (
        total == 16
    ), "Total usable sites should be {expected_total}, identified {x}".format(
        x=total, expected_total=expected_total
    )


def test_unique():
    test_df1 = DataFrame([["A"], ["B"]])
    test_df2 = DataFrame(["A"], ["A"])
    assert (
        is_unique(DataFrame(test_df1)) == False
    ), "Failed unique test, misidentified [A,B] as being composed of a single unique entity"
    assert (
        is_unique(DataFrame(test_df2)) == True
    ), "Failed unique test, misidentified [A,A] as being composed of a multiple unique entities"


def check_combinations(combo_barcodes):
    arr_2d = combo_barcodes.to_numpy()
    assert (
        np.array_equal(arr_2d[0], arr_2d[1]) == False
    ), "Duplicate barcodes were sampled"


def calculate_RH_sample(h_mono, h_poly):
    if h_mono > h_poly:
        rh = (h_mono - h_poly) / h_mono
    else:
        rh = (h_mono - h_poly) / h_mono
    return rh


def bootstrap(array, iterations):
    resample = np.random.choice(array, (iterations, len(array)), replace=True)
    return np.mean(resample, axis=1)


def rh_classifier(rh):
    rh = float(rh)
    if rh > 0.4:
        return "cotx"
    elif rh > 0.3:
        return "cotx_probable"
    elif rh > -0.1:
        return "coi=2"
    elif rh > -0.2:
        return "coi=2_probable"
    elif rh > -0.4:
        return "coi=3"
    elif rh > -0.6:
        return "coi=3_probable"
    elif rh > -0.8:
        return "coi=4_probable"
    else:
        return "coi=4"


def calculate_shannons_index(p_array):
    h = -np.sum([p_array * np.log(p_array)])
    return h


def calculate_evenness(h_array, k):
    e = h_array / np.log(k)
    return e


def calculate_h12(p_array):
    p_array = sorted(p_array)
    p1 = p_array[0]
    p2 = p_array[1]
    other_p = np.asarray(p_array[2:])
    sum_rest = np.sum(other_p**2)
    h12 = (p1 + p2) ** 2 + sum_rest
    return h12


def interpret_cotx_simbinbarcode(sim, cotx_event=1):
    """wokrks only for coi =2"""
    sim = sim[str(cotx_event)]
    converted_barcodes = []
    for binbarcode in sim:
        if binbarcode == 0:
            barcode = 24 * [0]
        else:
            barcode = [int(i) for i in list("{0:0b}".format(binbarcode))]
            while len(barcode) < 24:
                barcode = [
                    0
                ] + barcode  # the conversion does not preserve leading zeros
        converted_barcodes.append(barcode)
    return np.asarray(converted_barcodes)


def load_cotx_simulations():
    cotx_simulation_data = defaultdict(
        lambda: defaultdict(list)
    )  # dict[initial_coi][cotx_event]

    cotx_file = "cotx_barcodes_2.txt"
    cotx_data = json.load(open(cotx_file))
    for simulation in cotx_data:
        for n_cotx_event in [1, 2, 3]:
            if len(simulation[str(n_cotx_event)]) > 1:
                cotx_simulation_data[2]["cotx_event{n}".format(n=n_cotx_event)].append(
                    interpret_cotx_simbinbarcode(simulation, cotx_event=n_cotx_event)
                )

    for initial_coi in [3, 4, 5]:
        cotx_file = "cotx_barcodes_{x}.txt".format(x=initial_coi)
        cotx_data = json.load(open(cotx_file))
        for simulation in cotx_data:
            for n_cotx in [1, 2, 3]:
                if len(simulation[str(n_cotx)]) > 1:
                    cotx_simulation_data[initial_coi][
                        "cotx_event{x}".format(x=n_cotx)
                    ].append(np.asarray(simulation[str(n_cotx)]))
    return cotx_simulation_data


def wilson(p, n, z=1.96):
    denominator = 1 + z**2 / n
    centre_adjusted_probability = p + z * z / (2 * n)
    adjusted_standard_deviation = np.sqrt((p * (1 - p) + z * z / (4 * n)) / n)

    lower_bound = (
        centre_adjusted_probability - z * adjusted_standard_deviation
    ) / denominator
    upper_bound = (
        centre_adjusted_probability + z * adjusted_standard_deviation
    ) / denominator
    return lower_bound, upper_bound


class BarcodeStats:
    cpalette_converter = {
        "green": "viridis",
        "blue": "mako",
        "cornflowerblue": "mako",
        "crimson": "rocket",
        "orange": "Oranges",
        "purple": "Purples_r",
    }

    # Set up field names from sheet:
    multi_poly_field = "M_P"
    sample_name_field = "Sample_Name"

    def __init__(self, input_file, ISO3, barcode_file_path, sheet_name=None, adjusted_n=False):

        self.ISO3 = ISO3
        self.barcode_file_path = barcode_file_path

        # Ingest the input file:
        if input_file.endswith(".xls") or input_file.endswith(".xlsx"):
            if not sheet_name:
                self.master_df = DataFrame(read_excel(input_file, sheet_name=ISO3))
            else:
                self.master_df = DataFrame(read_excel(input_file, sheet_name=sheet_name))
        else:
            sep = "\t" if input_file.endswith(".tsv") else ","
            self.master_df = DataFrame(read_csv(input_file, sep=sep, header=0))

        # Ingest the barcode file:
        tmp_barcode_def_df = DataFrame(read_csv(barcode_file_path, sep="\t"))
        self.loci_position_names = list(tmp_barcode_def_df["name"].values)
        self.barcode_pos_dict = {name: f"{contig}:{pos}" for name, contig, pos in
                                 zip(tmp_barcode_def_df['name'], tmp_barcode_def_df['chr'], tmp_barcode_def_df['pos'])}
        # We don't have to, but we should clean up our mess:
        del tmp_barcode_def_df

        # Define all the fields we're going to create:
        self.poly_het_dict = None
        self.total_mono = None
        self.repeat_haps = None
        self.unique_mono_count = None
        self.popgen_stats = None
        self.haplotype_counts = None
        self.loci_allele_dict = None
        self.n_total = None
        self.n_poly = None
        self.n_singles = None
        self.poly_barcode_year_dfs = None
        self.mono_barcode_year_dfs = None
        self.barcode_year_dfs = None
        self.chrono_years = None
        self.mono_barcode_df = None
        self.barcodes_df = None
        self.poly_het_timeseries = None
        self.poly_barcode_het_dist = None
        self.poly_barcode_het_avg = None
        self.poly_samples = None
        self.RH_barcode_dict = None
        self.observed_RH = None
        self.H_mono_barcodes = None
        self.RH_df = None
        self.poly_df = None
        self.RH_yearly_averages = None
        self.RH_yearly_variances = None
        self.RH_average = None
        self.RH_weighted_var = None
        self.RH_ci = None
        self.model_expectations = None
        self.cotx_simulation_data = None

        # Now process our data:
        self.extract_region(ISO3)
        self.extract_mono_poly_stats()
        self.calc_stats()
        self.calculate_polyhet_timeseries()
        self.quantify_observed_poly_het(adjusted_n)
        self.calculate_RH()

    def label_haplotypes(self):
        barcode_haplotypes_dict = {}
        unique_barcode_haplotype = 0
        barcode_haplotypes = []
        for row in self.barcodes_df[self.barcodes_df.columns[7:31]].to_numpy():
            barcode = "".join(row)
            if barcode not in barcode_haplotypes_dict:
                barcode_haplotypes_dict[barcode] = unique_barcode_haplotype
                unique_barcode_haplotype += 1
            barcode_haplotypes.append(barcode_haplotypes_dict[barcode])
        self.barcodes_df["Haplotype"] = barcode_haplotypes

    def extract_region(self, ISO3):
        tmp_df = self.master_df[self.master_df["ISO3"] == ISO3]
        for position in self.loci_position_names:
            tmp_df[position] = [
                x.strip().upper() for x in tmp_df[position]
            ]  # noticed a tab formatting in these columns

        tmp_df[BarcodeStats.multi_poly_field] = [x.strip() for x in tmp_df[BarcodeStats.multi_poly_field]]

        control_samples = [
            "3D7",
            "3D7-1",
            "3D7-2",
            "3d7",
            "Dd2-MOD",
            "Dd2-Mod",
            "Dd2/Mod",
            "Dd2_MOD",
            "DD2",
            "Dd2",
        ]

        self.barcodes_df = tmp_df[
            (tmp_df["X"] <= 2) & (~tmp_df[BarcodeStats.sample_name_field].isin(control_samples))
        ]
        self.mono_barcode_df = self.barcodes_df[
            (self.barcodes_df["X"] <= 2)
            & (self.barcodes_df[BarcodeStats.multi_poly_field] == "M")
            & (~self.barcodes_df[BarcodeStats.sample_name_field].isin(control_samples))
        ]
        self.chrono_years = sorted(Counter(self.barcodes_df["Year"]).keys())
        self.label_haplotypes()

    def extract_mono_poly_stats(self, fail_threshold=2):
        self.barcode_year_dfs = {}
        self.mono_barcode_year_dfs = {}
        self.poly_barcode_year_dfs = {}
        for year in self.chrono_years:
            self.barcode_year_dfs[year] = self.barcodes_df[
                (self.barcodes_df["Year"] == year)
                & (self.barcodes_df["X"] <= fail_threshold)
            ]

            self.mono_barcode_year_dfs[year] = self.barcodes_df[
                (self.barcodes_df["Year"] == year)
                & (self.barcodes_df["X"] <= fail_threshold)
                & (self.barcodes_df[BarcodeStats.multi_poly_field] == "M")
            ]

            self.poly_barcode_year_dfs[year] = self.barcodes_df[
                (self.barcodes_df["Year"] == year)
                & (self.barcodes_df["X"] <= fail_threshold)
                & (self.barcodes_df[BarcodeStats.multi_poly_field] == "P")
            ]
            self.mono_barcode_year_dfs[year].reset_index(drop=True, inplace=True)
            self.poly_barcode_year_dfs[year].reset_index(drop=True, inplace=True)
        self.n_singles = np.asarray(
            [len(self.mono_barcode_year_dfs[year]) for year in self.chrono_years]
        )
        self.n_poly = np.asarray(
            [len(self.poly_barcode_year_dfs[year]) for year in self.chrono_years]
        )

        self.n_total = np.asarray(
            [len(self.barcode_year_dfs[year]) for year in self.chrono_years]
        )

        self.loci_allele_dict = {}
        for column in self.mono_barcode_df.columns[7:31]:
            counts = Counter(self.mono_barcode_df[column].to_numpy())
            counts.pop("X", 0)
            counts.pop("N", 0)
            major_allele = max(counts, key=counts.get)
            counts.pop(major_allele, 0)
            try:
                minor_allele = max(counts, key=counts.get)
            except:
                minor_allele = None
            self.loci_allele_dict[column] = (major_allele, minor_allele)

        self.haplotype_counts = {}
        for year in self.chrono_years:
            counts = Counter(
                self.barcode_year_dfs[year][self.barcode_year_dfs[year][BarcodeStats.multi_poly_field] == "M"][
                    "Haplotype"
                ]
            )
            self.haplotype_counts[year] = counts

    def calc_stats(self):
        self.popgen_stats = defaultdict(list)
        n = self.n_singles + self.n_poly
        p = self.n_poly / n
        q = 1.0 - p
        variances = p * q / n
        wilson_interval = [wilson(prop, n_samples) for prop, n_samples in zip(p, n)]

        self.popgen_stats["p_poly_fract"] = list(p)
        self.popgen_stats["var_poly_fract"] = list(variances)
        weighted_mean, weighted_var, weighted_ci = inverse_var_weight(p, p * q / n)
        self.popgen_stats["poly_fract_inv_var"] = (
            weighted_mean,
            weighted_var,
            weighted_ci,
        )
        self.popgen_stats["poly_wilson"] = wilson_interval
        x = np.asarray(self.chrono_years)
        X = sm.add_constant(x)
        y = np.array(p)
        model = sm.OLS(y, X).fit()
        self.popgen_stats["poly_fract_model"] = model

        # print('Unique Mono Fract')
        self.unique_mono_count = {}
        self.repeat_haps = defaultdict(lambda: defaultdict(lambda: 0))
        for year in self.chrono_years:
            sorted_hid = sorted(
                self.haplotype_counts[year],
                key=lambda x: self.haplotype_counts[year][x],
                reverse=True,
            )
            for hid in sorted_hid:
                if self.haplotype_counts[year][hid] != 1:
                    self.repeat_haps[year][hid] = self.haplotype_counts[year][hid]
                else:
                    self.repeat_haps[year]["unique"] += 1
            self.unique_mono_count[year] = self.repeat_haps[year]["unique"]

        self.total_mono = np.asarray(
            [
                np.sum(list(self.haplotype_counts[year].values()))
                for year in self.chrono_years
            ]
        )
        p_mono_unique = (
            np.asarray([self.unique_mono_count[year] for year in self.chrono_years])
            / self.total_mono
        )
        p_mono_clonal = 1.0 - p_mono_unique
        var_mono_unique = (p_mono_unique * (1.0 - p_mono_unique)) / self.total_mono
        self.popgen_stats["p_mono_unique"] = list(p_mono_unique)
        self.popgen_stats["var_mono_unique"] = list(var_mono_unique)
        self.popgen_stats["wilson_mono_unique"] = [
            wilson(p, n) for p, n in zip(p_mono_unique, self.total_mono)
        ]
        x = np.asarray(self.chrono_years)
        X = sm.add_constant(x)
        y = np.array(p_mono_unique)
        model = sm.OLS(y, X).fit()
        self.popgen_stats["mono_unique_model"] = model
        weighted_mean, weighted_var, weighted_ci = inverse_var_weight(
            p_mono_unique, var_mono_unique
        )
        self.popgen_stats["mono_unique_inv_var"] = (
            weighted_mean,
            weighted_var,
            weighted_ci,
        )

        self.popgen_stats["p_mono_clonal"] = list(p_mono_clonal)
        self.popgen_stats["var_mono_clonal"] = list(
            var_mono_unique
        )  # it's the same because it is the inverse
        self.popgen_stats["wilson_mono_clonal"] = [
            wilson(p, n) for p, n in zip(p_mono_clonal, self.total_mono)
        ]
        x = np.asarray(self.chrono_years)
        X = sm.add_constant(x)
        y = np.array(p_mono_clonal)
        model = sm.OLS(y, X).fit()
        self.popgen_stats["mono_clonal_model"] = model
        weighted_mean, weighted_var, weighted_ci = inverse_var_weight(
            p_mono_clonal, var_mono_unique
        )
        self.popgen_stats["mono_clonal_inv_var"] = (
            weighted_mean,
            weighted_var,
            weighted_ci,
        )

        # calculate mono diversity - all
        for year in self.chrono_years:
            hap_ids = np.asarray(list(self.haplotype_counts[year].keys()))
            hap_counts = np.asarray(list(self.haplotype_counts[year].values()))
            hap_freqs = hap_counts / np.sum(hap_counts)
            sampling_idxes = np.random.choice(hap_ids, p=hap_freqs, size=(200, 200))
            shannon_idxes, evenness_scores, H12_scores = [], [], []
            for sampling_idx in sampling_idxes:
                sampled_counts = Counter(sampling_idx)
                sampled_freqs = np.asarray(list(sampled_counts.values())) / 200

                H12 = calculate_h12(sampled_freqs)
                shannon_idx = calculate_shannons_index(sampled_freqs)
                shannon_idxes.append(shannon_idx)
                evenness = calculate_evenness(
                    shannon_idx, len(list(sampled_counts.keys()))
                )
                evenness_scores.append(evenness)
                H12_scores.append(H12)

            self.popgen_stats["shannon_idx_mean"].append(np.mean(shannon_idxes))
            self.popgen_stats["evenness_mean"].append(np.mean(evenness_scores))
            self.popgen_stats["H12_mean"].append(np.mean(H12_scores))

            self.popgen_stats["shannon_idx_var"].append(np.var(shannon_idxes))
            self.popgen_stats["evenness_var"].append(np.var(evenness_scores))
            self.popgen_stats["H12_var"].append(np.var(H12_scores))

        self.popgen_stats["shannon_idx_mean"] = np.array(
            self.popgen_stats["shannon_idx_mean"]
        )
        model = sm.OLS(self.popgen_stats["shannon_idx_mean"], X).fit()
        self.popgen_stats["shannon_idx_model"] = model

        self.popgen_stats["evenness_mean"] = np.array(
            self.popgen_stats["evenness_mean"]
        )
        model = sm.OLS(self.popgen_stats["evenness_mean"], X).fit()
        self.popgen_stats["evenness_model"] = model

        self.popgen_stats["H12_mean"] = np.array(self.popgen_stats["H12_mean"])
        model = sm.OLS(self.popgen_stats["H12_mean"], X).fit()
        self.popgen_stats["H12_model"] = model

        if "mccoil_median" in self.barcodes_df.columns:
            for year, df in self.barcodes_df.groupby("Year"):
                self.popgen_stats["mccoil_coi"].append(np.mean(df["mccoil_median"]))
                self.popgen_stats["mccoil_coi_std"].append(np.std(df["mccoil_median"]))
            for year, df in self.barcodes_df.groupby("Year"):
                self.popgen_stats["mccoil_coi_poly"].append(
                    np.mean([x for x in df["mccoil_median"] if x >= 2])
                )
                self.popgen_stats["mccoil_coi_poly_std"].append(
                    np.std([x for x in df["mccoil_median"] if x >= 2])
                )

    def calculate_polyhet_timeseries(self):
        self.poly_het_dict = defaultdict(dict)
        for year in self.chrono_years:
            for position in self.loci_position_names:
                counts = Counter(self.poly_barcode_year_dfs[year][position].to_list())
                n_missing = counts.pop("X", 0)
                total = np.sum(list(counts.values()))
                n_het = counts.pop("N", 0)
                p_het = n_het / total
                # print(year, position, p_het, n_missing, n_missing/total, total)
                self.poly_het_dict[year][position] = (p_het, total)

        self.poly_het_timeseries = defaultdict(list)
        for loci in self.loci_position_names:
            for year in self.chrono_years:
                self.poly_het_timeseries[loci].append(self.poly_het_dict[year][loci][0])

    def quantify_observed_poly_het(self, adjustedN=False):
        self.poly_barcode_het_dist = defaultdict(list)
        self.poly_barcode_het_avg = {}
        self.poly_samples = {}
        if not adjustedN:
            # counting N from barcode
            for year in self.chrono_years:
                for i, row in enumerate(
                    self.poly_barcode_year_dfs[year][
                        self.poly_barcode_year_dfs[year].columns[7:31]
                    ].to_numpy()
                ):
                    barcode = "".join(row)
                    barcode_counts = Counter(barcode)
                    barcode_counts.pop("X", 0)
                    total = np.sum(list(barcode_counts.values()))
                    het = barcode_counts.pop("N")
                    H_poly_barcode = het / total
                    self.poly_barcode_het_dist[year].append(H_poly_barcode)
                self.poly_samples[year] = list(
                    self.poly_barcode_year_dfs[year][BarcodeStats.sample_name_field]
                )
                self.poly_barcode_het_avg[year] = np.mean(
                    self.poly_barcode_het_dist[year]
                )
        else:
            for year in self.chrono_years:
                total = 24.0 - np.asarray(self.poly_barcode_year_dfs[year]["X"])
                het = (
                    np.asarray(self.poly_barcode_year_dfs[year]["Adjusted_Het"]) / total
                )
                self.poly_barcode_het_dist[year] = het
                self.poly_samples[year] = list(
                    self.poly_barcode_year_dfs[year][BarcodeStats.sample_name_field]
                )
                self.poly_barcode_het_avg[year] = np.mean(het)

    def sample_poly_barcodes(self, year, coi=2, samples=100):
        combinations = np.random.randint(
            self.mono_barcode_year_dfs[year].index[0],
            self.mono_barcode_year_dfs[year].index[-1],
            size=(samples, coi),
        )
        sampled_poly_barcodes = []
        for i, combo in enumerate(combinations):
            sampled_haplotypes = self.mono_barcode_year_dfs[year].loc[
                combo, ["Haplotype"]
            ]
            flag = len(np.unique(list(sampled_haplotypes["Haplotype"]))) != coi
            while flag:
                combo = np.random.randint(
                    self.mono_barcode_year_dfs[year].index[0],
                    self.mono_barcode_year_dfs[year].index[-1],
                    coi,
                )
                combinations[i] = combo
                sampled_haplotypes = self.mono_barcode_year_dfs[year].loc[
                    combo, ["Haplotype"]
                ]
                flag = len(np.unique(list(sampled_haplotypes["Haplotype"]))) != coi
            combo_barcodes = self.mono_barcode_year_dfs[year].loc[
                combo, self.loci_position_names
            ]
            check_combinations(combo_barcodes)
            sampled_poly_barcodes.append(
                combo_barcodes.apply(lambda x: return_stats(x), axis=0).to_list()
            )

        sampled_poly_barcodes = np.asarray(sampled_poly_barcodes)

        return sampled_poly_barcodes

    def simulator(self, n_poly, n_iterations=1000, f_true=0, axis=0, coi=2):
        run_fn_tests()
        poly_simulations = defaultdict(list)
        for year, n in zip(self.chrono_years, n_poly):
            for n_rep in range(n_iterations):
                attempt_count = 1
                if n_rep % 100 == 0:
                    print(year, n, n_rep)
                b = self.sample_poly_barcodes(year, coi=coi, samples=n)
                heterozygotes = quantify_het(b, f_true, axis=axis)
                total = quantify_total(b, axis=axis)
                p_het = heterozygotes / total
                while (
                    np.isfinite(p_het).all() == False
                ):  # if zero is found in the total
                    assert attempt_count <= 3, "maximum number of attempts reached"
                    print("attempting resample {x}".format(x=attempt_count))
                    b = self.sample_poly_barcodes(year, samples=n)
                    heterozygotes = quantify_het(b, f_true, axis=axis)
                    total = quantify_total(b, axis=axis)
                    p_het = heterozygotes / total
                    attempt_count += 1
                poly_simulations[year].append(p_het)
            poly_simulations[year] = np.asarray(poly_simulations[year])
        return poly_simulations

    def sample_cotx_barcodes_from_mono(
        self, year, coi=2, samples=100, initial_coi=2, cotx_event=1
    ):
        # cotx_sim stats
        cotx_sim_filtered_data = [
            x
            for x in self.cotx_simulation_data[initial_coi][
                "cotx_event{x}".format(x=cotx_event)
            ]
            if len(x) >= coi
        ]
        random_sim_idxes = np.asarray(
            random.sample(range(len(cotx_sim_filtered_data)), samples)
        )
        sampled_cotx_barcodes = np.asarray(cotx_sim_filtered_data, dtype="object")[
            random_sim_idxes
        ]

        # superinfection layer
        combinations = np.random.randint(
            self.mono_barcode_year_dfs[year].index[0],
            self.mono_barcode_year_dfs[year].index[-1],
            size=(samples, initial_coi),
        )
        sampled_poly_barcodes = []
        for i, combo in enumerate(combinations):
            sampled_haplotypes = self.mono_barcode_year_dfs[year].loc[
                combo, ["Haplotype"]
            ]
            flag = len(np.unique(list(sampled_haplotypes["Haplotype"]))) != initial_coi
            while flag:
                combo = np.random.randint(
                    self.mono_barcode_year_dfs[year].index[0],
                    self.mono_barcode_year_dfs[year].index[-1],
                    initial_coi,
                )
                combinations[i] = combo
                sampled_haplotypes = self.mono_barcode_year_dfs[year].loc[
                    combo, ["Haplotype"]
                ]
                flag = (
                    len(np.unique(list(sampled_haplotypes["Haplotype"]))) != initial_coi
                )
            combo_barcodes = self.mono_barcode_year_dfs[year].loc[
                combo, self.loci_position_names
            ]
            check_combinations(combo_barcodes)

            relatedness_maps = []
            for n in range(coi):
                relatedness_maps.append(sampled_cotx_barcodes[i][n])

            combo_barcodes = combo_barcodes.to_numpy()
            cotx_strains = []
            for relatedness_map in relatedness_maps:
                tmp = [
                    combo_barcodes[strain_choice][position]
                    for position, strain_choice in enumerate(relatedness_map)
                ]
                cotx_strains.append(tmp)

            df = DataFrame(cotx_strains)  # [cotx_strain1, cotx_strain2])
            sampled_poly_barcodes.append(
                df.apply(lambda x: return_stats(x), axis=0).to_list()
            )

        sampled_poly_barcodes = np.asarray(sampled_poly_barcodes)

        return sampled_poly_barcodes

    def simulator_cotx(
        self, n_poly, n_iterations=1000, axis=0, coi=2, initial_coi=2, cotx_event=1
    ):
        run_fn_tests()
        poly_simulations = defaultdict(list)
        for year, n in zip(self.chrono_years, n_poly):
            for n_rep in range(n_iterations):
                attempt_count = 1
                if n_rep % 100 == 0:
                    print(year, n, n_rep)
                b = self.sample_cotx_barcodes_from_mono(
                    year,
                    coi=coi,
                    samples=n,
                    initial_coi=initial_coi,
                    cotx_event=cotx_event,
                )
                heterozygotes = quantify_het(b, Ftrue=0, axis=axis)
                total = quantify_total(b, axis=axis)
                p_het = heterozygotes / total
                while (
                    np.isfinite(p_het).all() == False
                ):  # if zero is found in the total
                    assert attempt_count <= 3, "maximum number of attempts reached"
                    print("attempting resample {x}".format(x=attempt_count))
                    b = self.sample_poly_barcodes(year, samples=n)
                    heterozygotes = quantify_het(b, Ftrue, axis=axis)
                    total = quantify_total(b, axis=axis)
                    p_het = heterozygotes / total
                    attempt_count += 1
                poly_simulations[year].append(p_het)
            poly_simulations[year] = np.asarray(poly_simulations[year])
        return poly_simulations

    def calculate_RH(self, n_poly_per_year=20, n_iter=20):
        self.RH_barcode_dict = {}
        self.observed_RH = defaultdict(list)

        print("Simulating mono barcode sampling for RH")
        H_mono_barcodes = self.simulator(
            [n_poly_per_year for x in self.n_poly], n_iter, 0, axis=1
        )
        self.calculate_RH_barcode_distribution(H_mono_barcodes)
        self.calculate_RHyear_distribution(H_mono_barcodes)
        self.H_mono_barcodes = H_mono_barcodes

        # actual samples
        minimum_sample = np.min(self.n_poly)
        for year in self.chrono_years:
            for i, H_poly_barcode in enumerate(self.poly_barcode_het_dist[year]):
                RH_sample_dist = [
                    calculate_RH_sample(np.mean(sim_trace), H_poly_barcode)
                    for sim_trace in self.H_mono_barcodes[year]
                ]
                self.observed_RH[year].append(RH_sample_dist)

        X = sm.add_constant(self.chrono_years)
        y = np.array(
            [np.mean(self.RH_barcode_dict[year]) for year in self.chrono_years]
        )
        model1 = sm.OLS(y, X).fit()
        self.RH_barcode_dict["model"] = model1

        data, Rh_sample_averages, cotx_averages = (
            [],
            defaultdict(list),
            defaultdict(list),
        )
        for year in self.chrono_years:
            for sample_name, RH_sample_dist in zip(
                self.poly_samples[year], self.observed_RH[year]
            ):
                data.append([sample_name, year, np.mean(RH_sample_dist)])
        self.RH_df = DataFrame(data)
        self.RH_df.columns = ["Sample", "Year", "RH"]
        self.RH_df["classification"] = self.RH_df["RH"].apply(rh_classifier)
        self.poly_df = pd.merge(
            self.master_df, self.RH_df, left_on=BarcodeStats.sample_name_field, right_on="Sample"
        )

        for year, df in self.RH_df.groupby("Year"):
            total_sample = np.asarray(df["RH"])
            cotx_total_samples = np.asarray(df["classification"])

            sampling_idxes = np.random.randint(0, len(total_sample), size=(100, 200))
            for sampling_idx in sampling_idxes:
                RH_sample = total_sample[sampling_idx]
                Rh_sample_averages[year].append(np.mean(RH_sample))

                cotx_samples = cotx_total_samples[sampling_idx]
                cotx_counts = Counter(cotx_samples)
                p_cotx = (cotx_counts["cotx"] + cotx_counts["cotx_probable"]) / len(
                    sampling_idx
                )
                cotx_averages[year].append(p_cotx)
        averages = np.asarray(
            [np.mean(Rh_sample_averages[year]) for year in self.chrono_years]
        )
        variances = np.asarray(
            [np.var(Rh_sample_averages[year]) for year in self.chrono_years]
        )
        weighted_mean, weighted_var, weighted_ci = inverse_var_weight(
            averages, variances
        )

        self.RH_yearly_averages = averages
        self.RH_yearly_variances = variances
        self.RH_average = weighted_mean
        self.RH_weighted_var = weighted_var
        self.RH_ci = weighted_ci

    def calculate_RH_barcode_distribution(self, poly_simulations):
        def distance_RH_trace(RH, H_mono_trace, H_poly_trace):
            H_mono_trace = np.asarray(H_mono_trace)
            H_poly_trace = np.asarray(H_poly_trace)
            distance = np.sum((H_mono_trace - (H_mono_trace * RH) - H_poly_trace) ** 2)
            return distance

        H_poly_trace = list(self.poly_barcode_het_avg.values())
        n_reps = len(poly_simulations[list(poly_simulations.keys())[0]])
        barcode_het_timetraces = []
        for irep in range(n_reps):
            timetrace = []
            for year in self.chrono_years:
                timetrace.append(np.mean(poly_simulations[year][irep]))
            barcode_het_timetraces.append(timetrace)

        RH_barcode_distribution = []
        for itrace in range(n_reps):
            output = scipy.optimize.minimize(
                lambda RH: distance_RH_trace(
                    RH, barcode_het_timetraces[itrace], H_poly_trace
                ),
                x0=[0.5],
                bounds=[(-1, 1)],
            )
            RH_barcode_distribution.append(output)
        self.RH_barcode_dict["total"] = [
            output.x[0] for output in RH_barcode_distribution
        ]

    def calculate_RHyear_distribution(self, poly_simulations):
        def distance_RHyear_trace(
            RH, H_mono_traces_dict, itrace, year, H_poly_dict=self.poly_barcode_het_dist
        ):
            distance = 0
            H_mono_trace = np.mean(H_mono_traces_dict[year][itrace])
            H_poly_trace = np.mean(H_poly_dict[year])
            distance_array = (H_mono_trace - (H_mono_trace * RH) - H_poly_trace) ** 2
            return np.sum(distance_array)

        for year in self.chrono_years:
            n_reps = len(poly_simulations[list(poly_simulations.keys())[0]])
            RH_year_distribution = []
            for itrace in range(n_reps):
                output = scipy.optimize.minimize(
                    lambda RH: distance_RHyear_trace(
                        RH, poly_simulations, itrace, year
                    ),
                    x0=[0.5],
                    bounds=[(0, 1)],
                )
                RH_year_distribution.append(output)

            self.RH_barcode_dict[year] = [
                output.x[0] for output in RH_year_distribution
            ]

    def simulate_coi_cotx_sweep(self, n_poly=200, n_iter=200, oocyst_alpha=2.5):
        self.model_expectations = defaultdict(lambda: defaultdict(list))
        self.cotx_simulation_data = load_cotx_simulations(oocyst_alpha)

        H_mono_barcode_coi = defaultdict(dict)
        for coi in [2, 3, 4, 5]:
            print("Simulating COI={coi} expectation".format(coi=coi))
            H_mono_barcode_coi[coi] = self.simulator(
                [n_poly for x in self.n_poly], n_iter, 0, axis=1, coi=coi
            )

        for year in self.chrono_years:
            for coi in [2, 3, 4, 5]:
                for sim_iteration in H_mono_barcode_coi[coi][year]:
                    mean_het = np.mean(self.H_mono_barcodes[year])
                    for sim in sim_iteration:
                        self.model_expectations[
                            "coi={x}, oocyst_alpha={alpha}".format(
                                x=coi, alpha=oocyst_alpha
                            )
                        ][year].append(calculate_RH_sample(mean_het, sim))

        for initial_coi in [2, 3, 4, 5]:
            for cotx_round in [1, 2, 3]:
                print(
                    "Initial COI = {initial_coi}, Simulating cotx round {x}".format(
                        initial_coi=initial_coi, x=cotx_round
                    )
                )
                H_mono_barcode_cotx = self.simulator_cotx(
                    [n_poly for x in self.n_poly],
                    n_iter,
                    axis=1,
                    coi=2,
                    initial_coi=initial_coi,
                    cotx_event=cotx_round,
                )
                for year in self.chrono_years:
                    for sim_iteration in H_mono_barcode_cotx[year]:
                        mean_het = np.mean(self.H_mono_barcodes[year])
                        for sim in sim_iteration:
                            self.model_expectations[
                                "cotx_{initial_coi}_{alpha}_{cotx_round}".format(
                                    initial_coi=initial_coi,
                                    alpha=oocyst_alpha,
                                    cotx_round=cotx_round,
                                )
                            ][year].append(calculate_RH_sample(mean_het, sim))

    def plot_sample_distribution(
        self, color, ax=None, x_annotate=-0.2, y_annotate=0.1, legend=True, title=None
    ):
        if not ax:
            fig, ax = plt.subplots()
        ax.bar(
            [year for year in self.chrono_years],
            self.n_singles,
            color="grey",
            alpha=0.3,
        )

        bar = ax.bar(
            [year for year in self.chrono_years],
            self.n_poly,
            color=color,
            bottom=self.n_singles,
        )
        for x, y, z in zip(self.chrono_years, self.n_singles, self.n_singles):
            x = round(x, 2)
            y = round(y, 2)
            ax.annotate(
                str(z),
                (x + x_annotate, y + y_annotate),
                fontsize=15,
                color="black",
                fontweight="bold",
            )

        for x, y, z in zip(
            self.chrono_years, self.n_singles + self.n_poly, self.n_poly
        ):
            x = round(x, 2)
            y = round(y, 2)
            ax.annotate(
                str(z),
                (x + x_annotate, y + y_annotate),
                fontsize=12,
                color=color,
                fontweight="bold",
            )

        legend_elements = [
            Patch(facecolor="grey", edgecolor="black", label="Monogenomic"),
            Patch(facecolor=color, edgecolor="black", label="Polygenomic"),
        ]
        if legend:
            ax.legend(handles=legend_elements)
        if not title:
            ax.set_title("Sample Distribution", fontsize=20)
        else:
            ax.set_title(title, fontsize=20, loc="left")
        ax.set_xticks(self.chrono_years)
        ax.set_xticklabels(self.chrono_years)
        ax.tick_params(labelsize=15)
        return ax

    def plot_mono_poly_fract(
        self,
        color,
        ax=None,
        x_annotate=-0.02,
        y_annotate=0.02,
        annotate_color="black",
        title=None,
    ):
        if not ax:
            fig, ax = plt.subplots()
        p_mono = self.n_singles / self.n_total
        bar = ax.bar(
            [year for year in self.chrono_years], p_mono, color="grey", alpha=0.3
        )
        for b, z in zip(bar, [round(_, 2) for _ in p_mono]):
            x, y = b._x0, b._height
            # print(x,y)
            if str(x) == "nan":
                x = 0
            if str(y) == "nan":
                y = 0
            x = round(x, 2)
            y = round(y, 2)
            ax.annotate(
                str(z),
                (x + x_annotate, y + y_annotate),
                fontsize=12,
                color=annotate_color,
                fontweight="bold",
            )

        ax.bar(
            [year for year in self.chrono_years],
            self.n_poly / self.n_total,
            color=color,
            bottom=self.n_singles / self.n_total,
        )
        if not title:
            ax.set_title("Mono vs Poly Fraction", fontsize=20)
        else:
            ax.set_title(title, fontsize=20, loc="left")
        ax.tick_params(labelsize=15)
        ax.set_xlim(self.chrono_years[0] - 1, self.chrono_years[-1] + 1)
        ax.set_xticks(self.chrono_years)
        ax.set_xticklabels(self.chrono_years)

    def plot_mono_hap_sharing(
        self,
        color,
        ax=None,
        annotate_color="black",
        x_annotate=0.2,
        y_annotate=0.05,
        title=None,
    ):
        if not ax:
            fig, ax = plt.subplots()

        for year in self.chrono_years:
            bottom, unique = 0, 0
            sorted_hid = sorted(
                self.haplotype_counts[year],
                key=lambda x: self.haplotype_counts[year][x],
                reverse=True,
            )
            cpalette = sns.color_palette(
                BarcodeStats.cpalette_converter[color], len(self.repeat_haps[year])
            )
            total = np.sum(list(self.repeat_haps[year].values()))
            for i, hid in enumerate(self.repeat_haps[year]):
                if hid != "unique":
                    height = self.repeat_haps[year][hid] / total
                    ax.bar(
                        [year],
                        height,
                        bottom=bottom,
                        color=cpalette[i],
                        edgecolor="black",
                    )
                    bottom += height
            # shared_fracts.append(round(bottom,2))
            bar = ax.bar(
                [year],
                self.repeat_haps[year]["unique"] / total,
                bottom=bottom,
                color="grey",
                alpha=0.3,
            )

        shared_fracts = [
            round(_, 2) for _ in 1.0 - np.asarray(self.popgen_stats["p_mono_unique"])
        ]
        for x, y, z in zip(self.chrono_years, shared_fracts, shared_fracts):
            x = round(x, 2)
            y = round(y, 2)
            ax.annotate(
                str(z),
                (x - x_annotate, y + y_annotate),
                fontsize=12,
                color=annotate_color,
                fontweight="bold",
            )

        if not title:
            ax.set_title("Mono Clonality", fontsize=20)
        else:
            ax.set_title(title, fontsize=20, loc="left")
        ax.tick_params(labelsize=15)
        ax.set_xlim(self.chrono_years[0] - 1, self.chrono_years[-1] + 1)
        ax.set_xticks(self.chrono_years)
        ax.set_xticklabels(self.chrono_years)

    def plot_persistent_clones(
        self, color, ax=None, x_annotate=[-0.3, 0.2], y_annotate=0.1, title=None
    ):
        if not ax:
            fig, ax = plt.subplots()
        cpalette = sns.color_palette(BarcodeStats.cpalette_converter[color], 10)
        hap_stats = defaultdict(dict)
        total_clusters = 1
        for year in self.chrono_years:
            sorted_hid = sorted(
                self.haplotype_counts[year],
                key=lambda x: self.haplotype_counts[year][x],
                reverse=True,
            )
            for hid in sorted_hid:
                if self.haplotype_counts[year][hid] != 1:
                    hap_stats[hid][year] = self.haplotype_counts[year][hid]
                    total_clusters += 1
        cpalette = sns.color_palette(
            BarcodeStats.cpalette_converter[color], total_clusters
        )

        count = 1
        flipper = 0
        for haplotype in hap_stats:
            x_array, y_array, s_array = [], [], []
            for year in hap_stats[haplotype]:
                x_array.append(year)
                y_array.append(count)
                s_array.append(hap_stats[haplotype][year] * 100)
            ax.scatter(
                x_array, y_array, s_array, color=cpalette[count - 1], edgecolor="black"
            )
            ax.plot(x_array, y_array, color=cpalette[count - 1])
            for i, txt in enumerate(s_array):
                ax.annotate(
                    int(txt / 100.0),
                    (x_array[i] + x_annotate[flipper], y_array[i] - y_annotate),
                    fontsize=12,
                    color="black",
                    fontweight="bold",
                )
            if flipper == 0:
                flipper = 1
            else:
                flipper = 0
            count += 1
        if not title:
            ax.set_title("Mono Clonality", fontsize=20)
        else:
            ax.set_title(title, fontsize=20, loc="left")
        ax.tick_params(labelsize=15)
        ax.tick_params(left=False, labelleft=False)

        ax.set_xlim(self.chrono_years[0] - 1, self.chrono_years[-1] + 1)
        ax.set_xticks(self.chrono_years)
        ax.set_xticklabels(self.chrono_years)

    def plot_longitudinal(
        self, field, color="orange", ax=None, inverse_var=False, title=None
    ):
        if not ax:
            fig, ax = plt.subplots()
        fields = {
            "mono_unique": (
                "p_mono_unique",
                "var_mono_unique",
                "mono_unique_model",
                "mono_unique_inv_var",
            ),
            "poly_fract": (
                "p_poly_fract",
                "var_poly_fract",
                "poly_fract_model",
                "poly_fract_inv_var",
            ),
            "cotx": ("cotx_average", "cotx_var", "cotx_inv_var"),
        }

        p = self.popgen_stats[fields[field][0]]
        variances = self.popgen_stats[fields[field][1]]

        if not inverse_var:
            model = self.popgen_stats[fields[field][2]]
            x = np.asarray(self.chrono_years)
            X = sm.add_constant(x)
            ax.scatter(self.chrono_years, p, color=color)
            ax.plot(x, model.predict(X), color=color)
            ax.fill_between(
                self.chrono_years,
                p + 2.5 * np.sqrt(variances),
                p - 2.5 * np.sqrt(variances),
                alpha=0.3,
                color=color,
            )
            ax.set_ylim(0, 1)

            ax.set_title(field, fontsize=20)
            # ax.set_xlabel('Year', fontsize = 15)
            ax.set_ylabel("Proportion", fontsize=15)
            ax.tick_params(labelsize=15)
            ax.set_xlim(self.chrono_years[0] - 1, self.chrono_years[-1] + 1)
            ax.set_xticks(self.chrono_years)
            ax.set_xticklabels(self.chrono_years)

        else:
            shifted_x = [x + 1 for x in range(len(self.chrono_years))]
            ax.errorbar(
                shifted_x,
                p,
                color=color,
                yerr=1.96 * np.sqrt(variances),
                markersize=20,
                markeredgecolor=color,
                markerfacecolor=color,
                fmt=".",
                ecolor=color,
                capsize=10,
            )

            weighted_mean, weighted_var, weighted_coi = self.popgen_stats[
                fields[field][3]
            ]
            x = [shifted_x[0], shifted_x[-1]]
            ax.plot(x, [weighted_mean for _ in x], color=color, linewidth=3)
            ax.fill_between(
                x,
                [weighted_coi[0] for _ in x],
                [weighted_coi[1] for _ in x],
                color=color,
                linewidth=3,
                alpha=0.2,
            )
            ax.set_ylim(0, 1)
            if not title:
                ax.set_title(field, fontsize=20)
            else:
                ax.set_title(title, fontsize=20, loc="left")
            # ax.set_xlabel('Year', fontsize = 15)
            ax.set_ylabel("Proportion", fontsize=15)
            ax.tick_params(labelsize=15)
            ax.set_xlim(shifted_x[0] - 0.5, shifted_x[-1] + 0.5)
            ax.set_xticks(shifted_x)
            ax.set_xticklabels(shifted_x)

    def plot_RH_average_confidence(self, color="orange", ax=None):
        if not ax:
            fig, ax = plt.subplots()
        RH = np.mean(self.RH_barcode_dict["total"])
        RH_ci = (
            np.percentile(self.RH_barcode_dict["total"], 2.5),
            np.percentile(self.RH_barcode_dict["total"], 97.5),
        )
        ax.hist(self.RH_barcode_dict["total"], color="orange")

        ax.set_xlabel(r"$R_{H}$", fontsize=15)
        ax.set_ylabel("Freq", fontsize=15)

    def plot_RHsample_longitudinal_average(self, color="orange", ax=None):
        if not ax:
            fig, ax = plt.subplots()
        y = np.array(
            [np.mean(self.RH_barcode_dict[year]) for year in self.chrono_years]
        )
        X = sm.add_constant(self.chrono_years)

        ax.scatter(self.chrono_years, y, color=color)
        ax.plot(
            self.chrono_years, self.RH_barcode_dict["model"].predict(X), color=color
        )
        ax.boxplot(
            [self.RH_barcode_dict[year] for year in self.chrono_years],
            positions=self.chrono_years,
            showfliers=False,
            notch=True,
            patch_artist=True,
            boxprops=dict(facecolor=color, color=color),
            capprops=dict(color=color),
            whiskerprops=dict(color=color),
            flierprops=dict(color=color, markeredgecolor=color),
            medianprops=dict(color=color),
        )

        ax.set_ylim(0, 0.5)
        ax.tick_params(axis="both", labelsize=15)
        ax.set_ylabel(r"$R_{H}$", fontsize=20)
        ax.set_xlabel("Year", fontsize=20)

    def plot_RHsample_longitudinal(self, color="orange", ax=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(12, 5))
        b = sns.swarmplot(x="Year", y="RH", data=self.RH_df, color=color, ax=ax)
        ax.plot(
            [0 - 0.5, len(self.chrono_years) + 0.5],
            [self.RH_average, self.RH_average],
            color="black",
            linewidth=3,
        )
        ax.tick_params(axis="both", labelsize=15)

        sns.boxplot(
            showmeans=True,
            # meanline=True,
            meanprops={
                "marker": "s",
                "markerfacecolor": "white",
                "markeredgecolor": color,
                "markersize": "12",
            },
            medianprops={"visible": False},
            whiskerprops={"visible": False},
            zorder=10,
            x="Year",
            y="RH",
            data=self.RH_df,
            showfliers=False,
            showbox=False,
            showcaps=False,
            ax=ax,
        )
        ax.set_xlabel("Year", fontsize=20)
        ax.set_ylabel(r"$R_{H}$", fontsize=20)
        legend_elements = [
            Patch(
                facecolor="black",
                edgecolor="black",
                label=r"$R_{H}=$"
                + str(round(self.RH_average, 2))
                + " "
                + "({ci1},{ci2})".format(
                    ci1=str(round(self.RH_ci[0], 2)), ci2=str(round(self.RH_ci[1], 2))
                ),
            )
        ]

        ax.legend(handles=legend_elements, fontsize=15)
        ax.set_ylim(-1.1, 1.0)

    def plot_cotx_sweep(self, ax=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(12, 5))

        simulation_boxplot_results = []
        for year in self.chrono_years:
            for key in self.model_expectations:
                for RH in self.model_expectations[key][year]:
                    simulation_boxplot_results.append([year, RH, key])
        df = DataFrame(simulation_boxplot_results)
        df.columns = ["Year", "RH", "Condition"]
        df_melt = df.melt(id_vars=["Year", "Condition"], value_vars="RH")
        cotx_colors = sns.color_palette("rocket", 4)
        superinfection_colors = sns.color_palette("mako_r", 3)
        order = [
            "cotx_5_2_3",
            "cotx_4_2_3",
            "cotx_3_2_3",
            "cotx_2_2_3",
            "cotx_5_2_2",
            "cotx_4_2_2",
            "cotx_3_2_2",
            "cotx_2_2_2",
            "cotx_5_2_1",
            "cotx_4_2_1",
            "cotx_3_2_1",
            "cotx_2_2_1",
            "coi=2",
            "coi=3",
            "coi=4",
        ]
        colors = 3 * list(cotx_colors.as_hex()) + list(superinfection_colors.as_hex())
        custom_pal = {}
        for x, y in zip(order, colors):
            custom_pal[x] = y

        ax.fill_between([-0.5, 3.5], [1.0, 1.0], [-1.5, -1.5], color="grey", alpha=0.1)
        ax.fill_between(
            [8 - 0.5, 11.5], [1.0, 1.0], [-1.5, -1.5], color="grey", alpha=0.1
        )

        b = sns.boxplot(
            data=df_melt,
            x="Condition",
            y="value",
            showfliers=False,
            ax=ax,
            palette=custom_pal,
            showmeans=True,
            meanprops={
                "marker": "o",
                "markerfacecolor": "white",
                "markeredgecolor": "black",
                "markersize": "10",
            },
            order=order,
        )

        ax.tick_params(axis="both", labelsize=15)
        ax.set_ylabel(r"$R_{H}$", fontsize=20)
        ax.set_xlabel("Condition", fontsize=20)

        line3 = 4 * ["3"] + 4 * ["2"] + 4 * ["1"] + ["COI=2", "COI=3", "COI=4"]
        ax.set_xticklabels(line3, rotation=45)

        legend_elements = [
            Patch(
                facecolor=cotx_colors[0], edgecolor="black", label=r"$COI_{i,cotx}=5$"
            ),
            Patch(
                facecolor=cotx_colors[1], edgecolor="black", label=r"$COI_{i,cotx}=4$"
            ),
            Patch(
                facecolor=cotx_colors[2], edgecolor="black", label=r"$COI_{i,cotx}=3$"
            ),
            Patch(
                facecolor=cotx_colors[3], edgecolor="black", label=r"$COI_{i,cotx}=2$"
            ),
            Patch(
                facecolor=superinfection_colors[0],
                edgecolor="black",
                label=r"$COI_{i,super}=2$",
            ),
            Patch(
                facecolor=superinfection_colors[1],
                edgecolor="black",
                label=r"$COI_{i,super}=3$",
            ),
            Patch(
                facecolor=superinfection_colors[2],
                edgecolor="black",
                label=r"$COI_{i,super}=4$",
            ),
        ]

        ax.legend(handles=legend_elements, fontsize=12)

        ax.plot([0, 14], [0, 0], color="black", linestyle="--")
        ax.plot([0, 14], [0.3, 0.3], linestyle="--", color="crimson")
        ax.annotate("Cotx Detection\n     Threshold", [12.5, 0.35], color="crimson")

    def plot_RH_classification(self, color="orange", ax=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(12, 5))
        classification_counts = Counter(self.RH_df["classification"])
        total = np.sum(list(classification_counts.values()))
        x_array = np.asarray([0, 1, 2, 3])

        cat1 = [
            classification_counts[key] / total
            for key in ["cotx", "coi=2", "coi=3", "coi=4"]
        ]
        cat2 = [
            classification_counts[key] / total
            for key in [
                "cotx_probable",
                "coi=2_probable",
                "coi=3_probable",
                "coi=4_probable",
            ]
        ]
        ax.bar(x_array, cat1, color=color)
        ax.bar(x_array, cat2, bottom=cat1, color=color)

        proportions = np.asarray(cat1) + np.asarray(cat2)

        stdev_array = []
        for p in proportions:
            wilson_low, wilson_high = wilson(p, total)
            rel_low_boundary = p - wilson_low
            rel_high_boundary = wilson_high - p
            stdev_array.append((rel_low_boundary, rel_high_boundary))
        stdev_array = np.asarray(stdev_array).T

        ax.errorbar(
            x_array, proportions, yerr=stdev_array, fmt=".", capsize=5, color="black"
        )

        ax.set_xticks(x_array)
        ax.set_xticklabels(
            ["Cotransmission", "COI=2", "COI=3", "COI=4"], fontsize=15, rotation=45
        )
        ax.set_ylabel("Proportion", fontsize=15)
        ax.tick_params(labelsize=15)
        ax.set_ylim(0, 1)
        # ax.set_title('2020', fontsize = 20)

    def generate_summary_report(self, output_file=None):
        data_report = {}
        data_report["n_singles"] = self.n_singles
        data_report["n_poly"] = self.n_poly
        data_report["n_total"] = self.n_singles + self.n_poly

        data_report["poly_fract_data"] = [
            (round(p, 2), [round(x, 2) for x in ci])
            for p, ci in zip(
                self.popgen_stats["p_poly_fract"], self.popgen_stats["poly_wilson"]
            )
        ]
        data_report["p_mono_unique"] = [
            (round(p, 2), [round(x, 2) for x in ci])
            for p, ci in zip(
                self.popgen_stats["p_mono_unique"],
                self.popgen_stats["wilson_mono_unique"],
            )
        ]

        data_report["realmccoilcoi"] = [
            (round(p, 2), [max(round(p - 1.96 * std, 2), 0), round(p + 1.96 * std, 2)])
            for p, std in zip(
                self.popgen_stats["mccoil_coi"], self.popgen_stats["mccoil_coi_std"]
            )
        ]
        data_report["realmccoilcoi_poly"] = [
            (round(p, 2), [max(round(p - 1.96 * std, 2), 0), round(p + 1.96 * std, 2)])
            for p, std in zip(
                self.popgen_stats["mccoil_coi_poly"],
                self.popgen_stats["mccoil_coi_poly_std"],
            )
        ]

        data_report["RH_array"] = [
            (round(p, 2), [round(p - 1.96 * std, 2), round(p + 1.96 * std, 2)])
            for p, std in zip(
                self.RH_yearly_averages, np.sqrt(self.RH_yearly_variances)
            )
        ]
        data_report_df = pd.DataFrame.from_dict(data_report).T
        data_report_df.columns = self.chrono_years
        if output_file:
            data_report_df.to_csv(output_file)
        return data_report_df


def show_stats(ISO3, color):
    fig = plt.figure(figsize=(20, 15))
    axes = [fig.add_subplot(3, 3, i + 1) for i in range(0, 9)]

    BS[ISO3].plot_sample_distribution(color, axes[0])
    BS[ISO3].plot_mono_poly_fract(color, axes[1], x_annotate=0.0)
    BS[ISO3].plot_longitudinal("poly_fract", color, axes[2])
    axes[2].set_title("Poly Fraction", fontsize=15)

    BS[ISO3].plot_mono_hap_sharing(color, axes[3], x_annotate=0.35, y_annotate=0.03)
    BS[ISO3].plot_persistent_clones(color, ax=axes[4], x_annotate=[-0.5, 0.5])
    BS[ISO3].plot_longitudinal("mono_unique", color, axes[5])
    axes[5].set_title("Unique Mono Fraction", fontsize=15)

    BS[ISO3].plot_RHsample_longitudinal(ax=axes[6], color=color)
    BS[ISO3].plot_RH_classification(ax=axes[7], color=color)

    # Calculate some stats:
    counts = Counter(BS[ISO3].RH_df["classification"])
    total = np.sum(list(counts.values()))
    p = (counts["cotx"] + counts["cotx_probable"]) / total
    stdev = np.sqrt(p * (1 - p) / total)

    # We should save both a CSV and a TSV:
    base_name = ISO3.replace(":", ".")  # We have to replace ':' here for ease of use with filenames.
    fig.savefig(f"{base_name}_summary_figure.svg")
    fig.savefig(f"{base_name}_summary_figure.png")


if __name__ == "__main__":

    # Set up our CLI args:
    parser = argparse.ArgumentParser(
        description=f"Processes P. falciparum data from a spreadsheet into actionable information (e.g. CoI estimates)."
    )

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-f', '--input-file',
                               help='TSV/CSV/Excel file containing data to process',
                               required=True)
    requiredNamed.add_argument('-s', '--sheet',
                               help='Sheet name to process.',
                               required=True)
    requiredNamed.add_argument('-b', '--barcodes',
                               help='Barcode file to use.',
                               required=True)
    args = parser.parse_args()

    # Do some validation here:
    if (args.input_file.endswith(".xls") or args.input_file.endswith(".xlsx")) and not args.sheet:
        print("ERROR: You must provide a sheet name with an excel file input.", file=sys.stderr)
        sys.exit(1)

    # Do the work:
    BS = {}
    ISO3 = args.sheet.replace(".", ":")  # We have to replace '.' with ':' here because of Terra conventions.
    sheet_name = ISO3.replace(":", "_")

    BS[ISO3] = BarcodeStats(
        args.input_file, ISO3, args.barcodes, sheet_name=sheet_name, adjusted_n=False
    )

    show_stats(ISO3, color="crimson")

    base_name = ISO3.replace(":", ".")  # We have to replace ':' here for ease of use with filenames.
    BS[ISO3].generate_summary_report(f"{base_name}_summary.csv")
    BS[ISO3].mono_barcode_df.to_csv(f"{base_name}_mono_barcodes.csv")
    BS[ISO3].poly_df.to_csv(f"{base_name}_poly_barcodes.csv")
