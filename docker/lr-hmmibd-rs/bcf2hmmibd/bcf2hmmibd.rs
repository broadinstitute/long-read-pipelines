//! bcf2hmmibd — standalone BCF -> hmmIBD genotype-table converter.
//!
//! Emits the exact genotype matrix that `hmmibd-rs --from-bcf` feeds to its HMM,
//! serialised as an hmmIBD text table (header `chrom<TAB>pos<TAB>sample...`; each
//! data row `chrname<TAB>pos<TAB>call...`; `-1` = missing). This lets the pipeline
//! do the (expensive, correctness-critical) BCF read + filtering ONCE and pass a
//! small text artifact around, instead of re-reading the BCF for every IBD run.
//!
//! It is NOT a reimplementation: it calls the crate's own conversion functions
//! verbatim — `BcfFilterArgs::{new_from_toml_file,new_from_builtin}`,
//! `BcfGenotype::new_from_processing_bcf`, `Samples::from_args`, and
//! `BcfGenotype::into_genotype_siteinfo` — the identical code path that
//! `InputData::from_args` runs for `--from-bcf`. Therefore the table it writes,
//! when read back with `hmmibd-rs -i`, reproduces `--from-bcf` exactly.
//!
//! CLI: reuses hmmibd-rs's own `Arguments`, so the relevant flags are identical:
//!   -i/--data-file1 <BCF>         input BCF (required)
//!   -o/--output <FILE>            output table ("-" or omitted => stdout)
//!   --bcf-read-mode <MODE>        dominant-allele (default) | first-ploidy | each-ploidy
//!   --bcf-filter-config <TOML>    site/genotype filter config (default: built-in)
//!   --min-snp-sep <BP>            min inter-SNP spacing (default 5), applied as in --from-bcf
//!   -b/--bad-file <FILE>          sample ids to exclude (one per line)
//! Other hmmibd-rs flags are accepted and ignored (only the conversion runs).

use anyhow::Result;
use clap::Parser;
use hmmibd_rs::args::Arguments;
use hmmibd_rs::bcf::{BcfFilterArgs, BcfGenotype};
use hmmibd_rs::samples::Samples;
use std::io::{BufWriter, Write};

fn main() -> Result<()> {
    let args = Arguments::parse();

    // Build BCF filter args exactly as InputData::from_args does for --from-bcf:
    // an explicit --bcf-filter-config TOML if given, else the built-in defaults.
    let mut bcf_filter_args = match args.bcf_filter_config.as_ref() {
        Some(path) => BcfFilterArgs::new_from_toml_file(path)?,
        None => BcfFilterArgs::new_from_builtin()?,
    };
    // The CLI --absent-ad-as flag wins over the filter-config value (mirrors
    // InputData::from_args). Controls how a genotype with no FORMAT/AD AND no
    // FORMAT/GT is called (missing vs ref); absent-AD genotypes that DO have GT
    // are resolved from GT regardless of this.
    if let Some(absent_ad_as) = args.absent_ad_as {
        bcf_filter_args.absent_ad_as = absent_ad_as;
    }

    // Authoritative conversion: per-mode read (dominant/first/each-ploidy) +
    // iterative missingness filtering over nonmissing_rates + consolidation.
    let bcf_gt = BcfGenotype::new_from_processing_bcf(
        &args.bcf_read_mode,
        &bcf_filter_args,
        &args.data_file1,
    )?;

    // Same sample selection (honours --bad-file) as the IBD path.
    let valid_samples = Samples::from_args(&args, Some(&bcf_gt))?;
    let sample_names = valid_samples.v().clone();

    // Build the site-oriented genotype matrix + site info EXACTLY as
    // read_data_dominant_genotype does — this is the matrix --from-bcf feeds the
    // HMM (min_snp_sep spacing applied here, idempotent on a later -i re-read).
    let (geno, sites) = bcf_gt.into_genotype_siteinfo(&valid_samples, args.min_snp_sep)?;

    let nsites = geno.get_nrows();
    let nsam = geno.get_ncols();
    assert_eq!(nsam, sample_names.len(), "matrix cols != sample count");
    let chrname_vec = sites.chrname_vec();
    let chr_idx_vec = sites.chr_idx_vec();
    let chr_pos_vec = sites.chr_pos_vec();
    assert_eq!(nsites, chr_idx_vec.len(), "matrix rows != site count");

    let mut out: Box<dyn Write> = match args.output.as_deref() {
        Some(path) if path != "-" => Box::new(BufWriter::new(std::fs::File::create(path)?)),
        _ => Box::new(BufWriter::new(std::io::stdout())),
    };

    // header: chrom \t pos \t <sample names...>
    write!(out, "chrom\tpos")?;
    for s in &sample_names {
        write!(out, "\t{s}")?;
    }
    writeln!(out)?;

    // one row per site: chrname \t pos \t <calls; -1 = missing>
    use std::fmt::Write as _;
    let mut line = String::with_capacity(16 + nsam * 3);
    for i in 0..nsites {
        line.clear();
        let _ = write!(line, "{}\t{}", chrname_vec[chr_idx_vec[i]], chr_pos_vec[i]);
        for j in 0..nsam {
            let v = geno.get_at(i, j);
            if v == u8::MAX {
                line.push_str("\t-1");
            } else {
                let _ = write!(line, "\t{v}");
            }
        }
        line.push('\n');
        out.write_all(line.as_bytes())?;
    }
    out.flush()?;
    Ok(())
}
