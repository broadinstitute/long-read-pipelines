#!/usr/bin/env python3
"""Build a FLARE ref-panel file from gnomAD HGDP+1KG sample metadata."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import sys
from pathlib import Path

LOG = logging.getLogger(__name__)

# Common column names in gnomAD HGDP+1KG metadata releases.
SAMPLE_COLUMNS = ("sample", "Sample", "SAMPLE", "s", "SampleID", "sample_id")
POP_COLUMNS = (
    "population",
    "Population",
    "pop",
    "Pop",
    "population_name",
    "Population_Name",
    "super_population",
    "Super_Population",
)


def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz") or str(path).endswith(".bgz"):
        return gzip.open(path, "rt")
    return path.open()


def pick_column(fieldnames: list[str], candidates: tuple[str, ...]) -> str | None:
    lower_map = {f.lower(): f for f in fieldnames}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="gnomAD HGDP+1KG metadata table (TSV/CSV, optionally gzipped)",
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Output ref-panel path (whitespace-delimited, no header)",
    )
    parser.add_argument(
        "--sample-column",
        default=None,
        help="Metadata column for sample ID (auto-detected if omitted)",
    )
    parser.add_argument(
        "--population-column",
        default=None,
        help="Metadata column for population/panel label (auto-detected if omitted)",
    )
    parser.add_argument(
        "--delimiter",
        default=None,
        help="Field delimiter (default: tab, or sniffed from header)",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    with open_maybe_gzip(args.metadata) as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delimiter = args.delimiter or ("\t" if "\t" in sample.splitlines()[0] else ",")
        reader = csv.DictReader(fh, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError("Metadata file has no header row")

        sample_col = args.sample_column or pick_column(reader.fieldnames, SAMPLE_COLUMNS)
        pop_col = args.population_column or pick_column(reader.fieldnames, POP_COLUMNS)
        if not sample_col or not pop_col:
            raise ValueError(
                f"Could not detect sample/population columns in {reader.fieldnames}. "
                "Pass --sample-column and --population-column explicitly."
            )

        LOG.info("Using sample column '%s' and population column '%s'", sample_col, pop_col)
        rows: list[tuple[str, str]] = []
        seen: set[str] = set()
        for row in reader:
            sid = row[sample_col].strip()
            pop = row[pop_col].strip()
            if not sid or not pop:
                continue
            if sid in seen:
                continue
            seen.add(sid)
            rows.append((sid, pop))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as out:
        for sid, pop in rows:
            out.write(f"{sid}\t{pop}\n")

    LOG.info("Wrote %d samples to %s", len(rows), args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
