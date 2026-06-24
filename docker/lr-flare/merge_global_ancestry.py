#!/usr/bin/env python3
"""Merge per-chromosome FLARE global ancestry files into one SAIGE-ready TSV."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import sys
from collections import defaultdict
from pathlib import Path

LOG = logging.getLogger(__name__)


def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def read_global_anc(path: Path) -> tuple[list[str], dict[str, list[float]]]:
    with open_maybe_gzip(path) as fh:
        header_line = None
        for line in fh:
            if line.strip():
                header_line = line
                break
        if header_line is None:
            raise ValueError(f"Empty global ancestry file: {path}")

        header = header_line.strip().split()
        if header[0] != "SAMPLE":
            raise ValueError(f"Expected SAMPLE header in {path}, got: {header}")

        ancestries = header[1:]
        fractions: dict[str, list[float]] = defaultdict(list)

        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != len(header):
                raise ValueError(f"Malformed line in {path}: {line}")
            sample = parts[0]
            values = [float(x) for x in parts[1:]]
            fractions[sample].append(values)

    return ancestries, fractions


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        type=Path,
        help="Per-chromosome .global.anc.gz files",
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Output TSV path (sample x ancestry fractions)",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    if not args.inputs:
        raise ValueError("No input files provided")

    reference_ancestries: list[str] | None = None
    accumulated: dict[str, list[list[float]]] = defaultdict(list)

    for path in args.inputs:
        ancestries, by_sample = read_global_anc(path)
        if reference_ancestries is None:
            reference_ancestries = ancestries
        elif ancestries != reference_ancestries:
            raise ValueError(
                f"Ancestry columns differ in {path}: {ancestries} vs {reference_ancestries}"
            )
        for sample, value_lists in by_sample.items():
            accumulated[sample].extend(value_lists)

    assert reference_ancestries is not None

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(["SAMPLE", *reference_ancestries])
        for sample in sorted(accumulated):
            chr_values = accumulated[sample]
            if not chr_values:
                continue
            n_chr = len(chr_values)
            avg = [
                sum(vals[i] for vals in chr_values) / n_chr
                for i in range(len(reference_ancestries))
            ]
            writer.writerow([sample, *[f"{v:.6f}" for v in avg]])

    LOG.info(
        "Merged %d chromosome files for %d samples -> %s",
        len(args.inputs),
        len(accumulated),
        args.out,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
