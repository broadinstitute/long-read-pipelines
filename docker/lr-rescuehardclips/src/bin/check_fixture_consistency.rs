use std::path::Path;

use anyhow::{Context, anyhow, ensure};
use noodles::bam;
use noodles::sam::{
    alignment::{
        record::Cigar as _,
        record::QualityScores as _,
        record::cigar::op::Kind,
        RecordBuf,
    },
};

const FIXTURES: &[&str] = &[
    "tests/fixtures/test.dev.ONT.withHardClips.bam",
    "tests/fixtures/test.dev.ONT.withHardClips.usingGRCh38seq.bam",
];

#[derive(Default)]
struct FileSummary {
    total_records: usize,
    primary_records: usize,
    secondary_records: usize,
    supplementary_records: usize,
    hard_clipped_records: usize,
    soft_clipped_records: usize,
}

fn read_consuming_cigar_len(record: &RecordBuf) -> anyhow::Result<usize> {
    record.cigar().iter().try_fold(0usize, |len, op_res| {
        let op = op_res?;
        let consumes_read = matches!(
            op.kind(),
            Kind::Match
                | Kind::Insertion
                | Kind::SoftClip
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
        );

        if consumes_read {
            Ok(len + op.len())
        } else {
            Ok(len)
        }
    })
}

fn has_op_kind(record: &RecordBuf, kind: Kind) -> anyhow::Result<bool> {
    for op_res in record.cigar().iter() {
        if op_res?.kind() == kind {
            return Ok(true);
        }
    }

    Ok(false)
}

fn validate_record(record: &RecordBuf) -> anyhow::Result<()> {
    let read_name = record
        .name()
        .map(|name| String::from_utf8_lossy(name.as_ref()).into_owned())
        .unwrap_or_else(|| "<missing-name>".to_string());
    let seq_len = record.sequence().len();
    let qual_len = record.quality_scores().len();
    let cigar_len = read_consuming_cigar_len(record)?;

    ensure!(
        seq_len == qual_len,
        "record {} has inconsistent sequence and quality lengths: seq={}, qual={}",
        read_name,
        seq_len,
        qual_len,
    );
    ensure!(
        seq_len == cigar_len,
        "record {} has inconsistent sequence length and read-consuming CIGAR span: seq={}, cigar={}",
        read_name,
        seq_len,
        cigar_len,
    );

    Ok(())
}

fn summarize_bam(path: &Path) -> anyhow::Result<FileSummary> {
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(path)
        .with_context(|| format!("failed to open {}", path.display()))?;
    let header = reader.read_header()?;
    let mut summary = FileSummary::default();

    for result in reader.record_bufs(&header) {
        let record = result?;
        validate_record(&record)
            .with_context(|| format!("consistency check failed in {}", path.display()))?;

        summary.total_records += 1;
        if record.flags().is_secondary() {
            summary.secondary_records += 1;
        } else if record.flags().is_supplementary() {
            summary.supplementary_records += 1;
        } else {
            summary.primary_records += 1;
        }

        if has_op_kind(&record, Kind::HardClip)? {
            summary.hard_clipped_records += 1;
        }
        if has_op_kind(&record, Kind::SoftClip)? {
            summary.soft_clipped_records += 1;
        }
    }

    Ok(summary)
}

fn main() -> anyhow::Result<()> {
    for fixture in FIXTURES {
        let path = Path::new(fixture);
        let summary = summarize_bam(path)?;

        if summary.total_records == 0 {
            return Err(anyhow!("{} contains no records", path.display()));
        }

        println!(
            "{}\n  total={} primary={} secondary={} supplementary={} hardclip={} softclip={}",
            path.display(),
            summary.total_records,
            summary.primary_records,
            summary.secondary_records,
            summary.supplementary_records,
            summary.hard_clipped_records,
            summary.soft_clipped_records,
        );
    }

    Ok(())
}
