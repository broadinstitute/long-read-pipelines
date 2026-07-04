#!/usr/bin/env python3
"""Prepare study and reference VCFs for FLARE local ancestry inference."""

from __future__ import annotations

import argparse
import gzip
import logging
import subprocess
import sys
import tempfile
from pathlib import Path

LOG = logging.getLogger(__name__)


def run_cmd(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    LOG.info("Running: %s", " ".join(cmd))
    proc = subprocess.run(cmd, check=False, text=True, capture_output=True)
    if proc.stdout.strip():
        LOG.debug(proc.stdout.strip())
    if proc.stderr.strip():
        LOG.info(proc.stderr.strip())
    if check and proc.returncode != 0:
        msg = proc.stderr.strip() or proc.stdout.strip() or "unknown error"
        raise RuntimeError(
            f"Command failed ({proc.returncode}): {' '.join(cmd)}\n{msg}"
        ) from None
    return proc


def count_records(path: Path) -> int:
    """Stream-count variants without loading records into Python memory."""
    proc = subprocess.Popen(
        ["bcftools", "view", "-H", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert proc.stdout is not None
    count = 0
    for _ in proc.stdout:
        count += 1
    stderr = proc.stderr.read() if proc.stderr else ""
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(
            f"bcftools view failed counting records in {path}: {stderr.strip()}"
        )
    return count


def log_file_size(path: Path, label: str) -> None:
    if not path.exists():
        LOG.warning("%s: %s does not exist", label, path)
        return
    size_gb = path.stat().st_size / (1024**3)
    LOG.info("%s: %s (%.2f GiB)", label, path.name, size_gb)


def log_file_stats(path: Path, label: str) -> None:
    log_file_size(path, label)
    LOG.info("%s record count: %s", label, count_records(path))


def bcftools_index(path: Path) -> None:
    run_cmd(["bcftools", "index", "-c", "-f", str(path)])


def read_sex_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    sex_by_sample: dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"Invalid sex map line (expected 2 columns): {line}")
            sample, sex = parts[0], parts[1].strip().lower()
            sex_by_sample[sample] = sex
    return sex_by_sample


def require_fasta_index(ref_fasta: Path, ref_fasta_fai: Path) -> None:
    expected = Path(f"{ref_fasta}.fai")
    if expected.exists():
        return
    if not ref_fasta_fai.exists():
        raise FileNotFoundError(f"Missing FASTA index: {ref_fasta_fai}")
    expected.symlink_to(ref_fasta_fai)
    LOG.info("Linked %s -> %s", expected, ref_fasta_fai)


def normalize_vcf(vcf: Path, ref_fasta: Path, out_path: Path) -> None:
    run_cmd(
        [
            "bcftools",
            "norm",
            "-f",
            str(ref_fasta),
            "-m",
            "-any",
            "-Ob",
            "-o",
            str(out_path),
            str(vcf),
        ]
    )
    log_file_size(out_path, "normalized")


def filter_biallelic_snps(vcf: Path, out_path: Path) -> None:
    run_cmd(
        [
            "bcftools",
            "view",
            "-m2",
            "-M2",
            "-v",
            "snps",
            "-Ob",
            "-o",
            str(out_path),
            str(vcf),
        ]
    )


def filter_maf(vcf: Path, maf: float, out_path: Path) -> None:
    run_cmd(
        [
            "bcftools",
            "view",
            "-q",
            f"{maf}:minor",
            "-Ob",
            "-o",
            str(out_path),
            str(vcf),
        ]
    )


def intersect_vcfs(
    gt_vcf: Path,
    ref_vcf: Path,
    gt_out: Path,
    ref_out: Path,
    *,
    output_vcf: bool,
) -> None:
    bcftools_index(gt_vcf)
    bcftools_index(ref_vcf)
    out_args = ["-Oz"] if output_vcf else ["-Ob"]

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        run_cmd(
            [
                "bcftools",
                "isec",
                "-p",
                str(tmpdir),
                "-n=2",
                "-w1",
                str(gt_vcf),
                str(ref_vcf),
            ]
        )
        run_cmd(["bcftools", "view", *out_args, "-o", str(gt_out), str(tmpdir / "0000.vcf")])
        run_cmd(["bcftools", "view", *out_args, "-o", str(ref_out), str(tmpdir / "0001.vcf")])

    if not output_vcf:
        bcftools_index(gt_out)
        bcftools_index(ref_out)


def thin_vcf(input_vcf: Path, output_vcf: Path, min_spacing: int) -> None:
    open_fn = gzip.open if str(input_vcf).endswith(".gz") else open
    out_open = gzip.open if str(output_vcf).endswith(".gz") else open

    last_kept_pos: int | None = None
    last_chrom: str | None = None
    kept = 0

    with open_fn(input_vcf, "rt") as fin, out_open(output_vcf, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos_s = fields[0], fields[1]
            pos = int(pos_s)
            if chrom != last_chrom:
                last_chrom = chrom
                last_kept_pos = None
            if last_kept_pos is not None and (pos - last_kept_pos) < min_spacing:
                continue
            fout.write(line)
            last_kept_pos = pos
            kept += 1

    LOG.info("Thinned to %s records in %s", kept, output_vcf.name)


def fix_male_chrx_hets(
    vcf: Path,
    sex_by_sample: dict[str, str],
    chromosome: str,
    out_path: Path,
) -> None:
    if chromosome not in ("chrX", "X") or not sex_by_sample:
        run_cmd(["bcftools", "view", "-Ob", "-o", str(out_path), str(vcf)])
        return

    males = [s for s, sex in sex_by_sample.items() if sex in ("male", "m")]
    if not males:
        LOG.warning("is_chr_x set but sex map contains no male samples")
        run_cmd(["bcftools", "view", "-Ob", "-o", str(out_path), str(vcf)])
        return

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        males_file = tmpdir / "males.txt"
        males_file.write_text("\n".join(males) + "\n")

        non_par = tmpdir / "non_par.bcf"
        run_cmd(
            [
                "bcftools",
                "view",
                "-r",
                "chrX:2781480-155701382",
                "-Ob",
                "-o",
                str(non_par),
                str(vcf),
            ]
        )
        bcftools_index(non_par)

        fixed_non_par = tmpdir / "fixed_non_par.bcf"
        run_cmd(
            [
                "bcftools",
                "+setGT",
                str(non_par),
                "-Ob",
                "-o",
                str(fixed_non_par),
                "--",
                "-t",
                "het",
                "-n",
                "c:'1|1'",
                "-s",
                str(males_file),
            ]
        )
        bcftools_index(fixed_non_par)

        par_vcf = tmpdir / "par.bcf"
        run_cmd(
            [
                "bcftools",
                "view",
                "-r",
                "chrX:10001-2781479,chrX:155701383-156030895",
                "-Ob",
                "-o",
                str(par_vcf),
                str(vcf),
            ]
        )
        bcftools_index(par_vcf)

        run_cmd(
            [
                "bcftools",
                "concat",
                "-a",
                "-Ob",
                "-o",
                str(out_path),
                str(par_vcf),
                str(fixed_non_par),
            ]
        )


def validate_phased_gt(vcf: Path, max_records: int = 1000) -> None:
    proc = subprocess.Popen(
        ["bcftools", "view", "-H", str(vcf)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert proc.stdout is not None
    checked = 0
    try:
        for line in proc.stdout:
            fields = line.rstrip("\n").split("\t")
            for gt in fields[9:]:
                if not gt or gt == ".":
                    continue
                fmt = gt.split(":")[0]
                if "/" in fmt:
                    raise ValueError(f"Unphased genotype found (expected '|'): {fmt}")
                if "." in fmt:
                    raise ValueError(f"Missing allele in genotype: {fmt}")
            checked += 1
            if checked >= max_records:
                break
    finally:
        proc.stdout.close()
        if proc.poll() is None:
            proc.kill()
        stderr = proc.stderr.read() if proc.stderr else ""
        proc.wait()
        if proc.returncode not in (0, -9) and checked < max_records:
            raise RuntimeError(stderr.strip() or "bcftools view failed during GT validation")
    LOG.info("Validated phased GT format on %s records", checked)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gt", required=True, type=Path, help="Study BCF/VCF for one chromosome")
    parser.add_argument("--ref", required=True, type=Path, help="Reference VCF for one chromosome")
    parser.add_argument("--ref-fasta", required=True, type=Path, help="GRCh38 reference FASTA")
    parser.add_argument("--ref-fasta-fai", required=True, type=Path, help="FASTA index (.fai)")
    parser.add_argument("--chromosome", required=True, help="Chromosome name, e.g. chr20 or chrX")
    parser.add_argument("--out-prefix", required=True, type=Path, help="Output prefix for filtered VCFs")
    parser.add_argument("--maf", type=float, default=0.01, help="Minimum MAF in study cohort")
    parser.add_argument("--thin-bp", type=int, default=20_000, help="Minimum spacing between retained SNPs")
    parser.add_argument("--is-chr-x", action="store_true", help="Apply chrX-specific handling")
    parser.add_argument("--sample-sex-map", type=Path, default=None, help="Sample sex map (SAMPLE SEX)")
    parser.add_argument("--skip-validation", action="store_true", help="Skip phased GT validation")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    require_fasta_index(args.ref_fasta, args.ref_fasta_fai)
    sex_by_sample = read_sex_map(args.sample_sex_map)
    if args.is_chr_x and not sex_by_sample:
        LOG.warning("chrX run without sample_sex_map; skipping male het correction")

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    gt_out = Path(f"{args.out_prefix}.gt.flare.vcf.gz")
    ref_out = Path(f"{args.out_prefix}.ref.flare.vcf.gz")

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        gt_norm = tmpdir / "gt.norm.bcf"
        ref_norm = tmpdir / "ref.norm.bcf"
        gt_snps = tmpdir / "gt.snps.bcf"
        ref_snps = tmpdir / "ref.snps.bcf"
        gt_maf = tmpdir / "gt.maf.bcf"
        gt_isec = tmpdir / "gt.isec.vcf.gz"
        ref_isec = tmpdir / "ref.isec.vcf.gz"
        gt_thin = tmpdir / "gt.thin.vcf.gz"
        ref_thin = tmpdir / "ref.thin.vcf.gz"
        gt_chrx = tmpdir / "gt.chrx.bcf"

        normalize_vcf(args.gt, args.ref_fasta, gt_norm)
        normalize_vcf(args.ref, args.ref_fasta, ref_norm)

        if args.is_chr_x and sex_by_sample:
            fix_male_chrx_hets(gt_norm, sex_by_sample, args.chromosome, gt_chrx)
            gt_for_filter = gt_chrx
        else:
            gt_for_filter = gt_norm

        filter_biallelic_snps(gt_for_filter, gt_snps)
        filter_biallelic_snps(ref_norm, ref_snps)
        filter_maf(gt_snps, args.maf, gt_maf)
        log_file_stats(gt_maf, "study MAF-filtered")

        intersect_vcfs(gt_maf, ref_snps, gt_isec, ref_isec, output_vcf=True)
        thin_vcf(gt_isec, gt_thin, args.thin_bp)
        thin_vcf(ref_isec, ref_thin, args.thin_bp)

        # Re-intersect after thinning so site sets match exactly.
        intersect_vcfs(gt_thin, ref_thin, gt_out, ref_out, output_vcf=True)

    bcftools_index(gt_out)
    bcftools_index(ref_out)

    if not args.skip_validation:
        validate_phased_gt(gt_out)

    LOG.info("Wrote %s and %s", gt_out, ref_out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
