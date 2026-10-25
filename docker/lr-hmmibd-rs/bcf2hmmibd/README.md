# bcf2hmmibd — standalone BCF → hmmIBD genotype-table converter

`bcf2hmmibd` converts a BCF directly into the hmmIBD text genotype table, producing
the **exact** genotype matrix that `hmmibd-rs --from-bcf` feeds to its HMM. This
lets the pipeline do the (expensive, correctness-critical) BCF read + filtering
once and pass a small text artifact around, instead of re-reading the BCF for every
IBD run — and it replaces the hand-written awk conversion in the WDL, which had
drifted from `--from-bcf`'s real filtering.

## Why it cannot drift

It is **not a reimplementation**. It is a new binary added *in-tree* to the pinned
hmmibd-rs v0.1.7 source that calls the crate's own conversion functions verbatim:

```
BcfFilterArgs::{new_from_toml_file,new_from_builtin}   # same filter config
  -> BcfGenotype::new_from_processing_bcf              # read_dom + iterative
                                                        #   missingness + consolidate
  -> Samples::from_args                                # same sample selection (--bad-file)
  -> BcfGenotype::into_genotype_siteinfo               # same matrix + min_snp_sep spacing
```

This is the identical path `InputData::from_args` runs for `--from-bcf`
(`src/data.rs`). The table it writes, read back with `hmmibd-rs -i`, reproduces
`--from-bcf` byte-for-byte.

## Vendored files

- `bcf2hmmibd.rs` — the binary source (built as `src/bin/bcf2hmmibd.rs`).
- `0001-expose-SiteInfoRaw-getters.patch` — adds three one-line getters
  (`chrname_vec`/`chr_idx_vec`/`chr_pos_vec`) to `SiteInfoRaw` so the separate bin
  target can read the site info back out. The only change to upstream code.

The Dockerfile copies both in after checking out the pinned release, applies the
patch, and `cargo build --release` builds `bcf2hmmibd` alongside `hmmibd-rs`.

## Usage

```
bcf2hmmibd -i <input.bcf> -o <table.txt>
           [--bcf-read-mode dominant-allele|first-ploidy|each-ploidy]   # default dominant-allele
           [--bcf-filter-config <config.toml>]   # default: built-in filters
           [--min-snp-sep <bp>]                  # default 5 (same as --from-bcf)
           [-b <bad-samples.txt>]                # samples to exclude, one per line
```

Output: header `chrom<TAB>pos<TAB>sample...`, one data row per site
`chrname<TAB>pos<TAB>call...`, `-1` = missing. chrom is the BCF contig NAME
(hmmibd-rs reads string chrom; feed with `hmmibd-rs -i`, not the integer-only C
hmmIBD). (CLI reuses hmmibd-rs's own `Arguments`, so other hmmibd-rs flags are
accepted and ignored.)

NOTE: because the CLI is hmmibd-rs's own, `--bcf-filter-config` requires the
`--from-bcf` flag to also be present (a clap group constraint). `--from-bcf` is a
no-op for bcf2hmmibd — it always converts — so just pass it alongside a config:
`bcf2hmmibd --from-bcf -i in.bcf -o out.txt --bcf-filter-config cfg.toml`.

## Pipeline wiring

The `HmmIBD` WDL (`FilterVcfForHmmIBD` task, `output_format=hmmibd_table`) uses
bcf2hmmibd to convert each per-shard filtered BCF in-image (bcftools + bcf2hmmibd
co-located in `lr-hmmibd-rs`), so the big BCF never leaves the machine. It writes
a `BcfFilterArgs` TOML from the WDL knobs with `nonmissing_rates = []`: that
iterative site+SAMPLE missingness refinement needs whole-callset missingness,
which isn't available when the callset is sharded across input VCFs (a sample
dropped in one shard but not another would desync the per-shard columns and break
the stitch). All per-SITE filters stay on, so the per-shard table is byte-identical
to `--from-bcf` run on that shard under the same config. `StitchHmmIBDTables`
version-sorts chrom so contig names order correctly.

## Equivalence (verified)

Synthetic BCF exercising every filter (PASS/non-PASS, indel major/minor,
monomorphic, all-low-depth→missing, multiallelic, min_snp_sep close pair, two
contigs), 12 samples:

1. **bcf2hmmibd == `--from-bcf` (byte-exact).** `hmmibd-rs -i <bcf2hmmibd table>`
   produces md5-identical `hmm.txt` and `hmm_fract.txt` to `hmmibd-rs --from-bcf`
   on the same BCF.
2. **hmmibd-rs == C hmmIBD on the same table (numerically exact).** Feeding the
   table to the bundled C `hmmIBD` and to `hmmibd-rs -i` gives numerically
   identical output (max relative field diff 0; only cosmetic `log_p` exponent
   formatting differs).

## Notes / behavior to be aware of

- **Positions are 0-based** (bcf_reader `record.pos()`), matching `--from-bcf`'s
  internal coordinates. The old awk conversion emitted 1-based VCF `POS` — i.e. it
  was off-by-one from `--from-bcf`. Biologically negligible for IBD (min_snp_sep is
  relative; a 1 bp shift in genetic position is immaterial).
- **`major_minor_alleles_must_be_snps` is a no-op in hmmibd-rs v0.1.7** (a bug:
  `get_maf_and_sorted_allele_count` never stores the allele index, so the SNP check
  always compares allele 0 to itself and returns true). `bcf2hmmibd` reproduces
  this exactly (that is the point — exact equivalence to hmmibd-rs). To actually
  remove indel alleles, do it upstream with `bcftools norm -m-any` + a
  biallelic/SNP-only filter.
