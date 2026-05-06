use std::{
    fs,
    fs::File,
    io::{BufWriter, Write},
    num::NonZeroUsize,
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

use anyhow::{Context, Result, bail};
use clap::Parser;
use rust_htslib::bam::{
    self, CompressionLevel, Read,
    record::{Aux, AuxArray, Record},
};

const TARGET_TAG: &[u8] = b"mx";
const DEFAULT_PROGRESS_SECONDS: u64 = 60;

#[derive(Debug, Parser)]
#[command(
    name = "lr_drop_duptag",
    version,
    about = "Drop duplicate mx aux tags from BAM records"
)]
struct Cli {
    #[arg(long)]
    input: PathBuf,

    #[arg(long = "output-bam")]
    output_bam: PathBuf,

    #[arg(long)]
    mismatches: PathBuf,

    #[arg(long, default_value_t = NonZeroUsize::new(8).unwrap())]
    threads: NonZeroUsize,

    #[arg(long, default_value_t = DEFAULT_PROGRESS_SECONDS)]
    progress_seconds: u64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum OwnedAux {
    Char(u8),
    I8(i8),
    U8(u8),
    I16(i16),
    U16(u16),
    I32(i32),
    U32(u32),
    Float(u32),
    Double(u64),
    String(String),
    HexByteArray(String),
    ArrayI8(Vec<i8>),
    ArrayU8(Vec<u8>),
    ArrayI16(Vec<i16>),
    ArrayU16(Vec<u16>),
    ArrayI32(Vec<i32>),
    ArrayU32(Vec<u32>),
    ArrayFloat(Vec<u32>),
}

impl OwnedAux {
    fn from_aux(aux: Aux<'_>) -> Self {
        match aux {
            Aux::Char(v) => Self::Char(v),
            Aux::I8(v) => Self::I8(v),
            Aux::U8(v) => Self::U8(v),
            Aux::I16(v) => Self::I16(v),
            Aux::U16(v) => Self::U16(v),
            Aux::I32(v) => Self::I32(v),
            Aux::U32(v) => Self::U32(v),
            Aux::Float(v) => Self::Float(v.to_bits()),
            Aux::Double(v) => Self::Double(v.to_bits()),
            Aux::String(v) => Self::String(v.to_string()),
            Aux::HexByteArray(v) => Self::HexByteArray(v.to_string()),
            Aux::ArrayI8(values) => Self::ArrayI8(values.iter().collect()),
            Aux::ArrayU8(values) => Self::ArrayU8(values.iter().collect()),
            Aux::ArrayI16(values) => Self::ArrayI16(values.iter().collect()),
            Aux::ArrayU16(values) => Self::ArrayU16(values.iter().collect()),
            Aux::ArrayI32(values) => Self::ArrayI32(values.iter().collect()),
            Aux::ArrayU32(values) => Self::ArrayU32(values.iter().collect()),
            Aux::ArrayFloat(values) => Self::ArrayFloat(values.iter().map(f32::to_bits).collect()),
        }
    }

    fn push_to_record(&self, record: &mut Record, tag: &[u8]) -> Result<()> {
        match self {
            Self::Char(v) => record.push_aux(tag, Aux::Char(*v))?,
            Self::I8(v) => record.push_aux(tag, Aux::I8(*v))?,
            Self::U8(v) => record.push_aux(tag, Aux::U8(*v))?,
            Self::I16(v) => record.push_aux(tag, Aux::I16(*v))?,
            Self::U16(v) => record.push_aux(tag, Aux::U16(*v))?,
            Self::I32(v) => record.push_aux(tag, Aux::I32(*v))?,
            Self::U32(v) => record.push_aux(tag, Aux::U32(*v))?,
            Self::Float(bits) => record.push_aux(tag, Aux::Float(f32::from_bits(*bits)))?,
            Self::Double(bits) => record.push_aux(tag, Aux::Double(f64::from_bits(*bits)))?,
            Self::String(v) => record.push_aux(tag, Aux::String(v.as_str()))?,
            Self::HexByteArray(v) => record.push_aux(tag, Aux::HexByteArray(v.as_str()))?,
            Self::ArrayI8(values) => {
                let aux_array: AuxArray<'_, i8> = values.as_slice().into();
                record.push_aux(tag, Aux::ArrayI8(aux_array))?;
            }
            Self::ArrayU8(values) => {
                let aux_array: AuxArray<'_, u8> = values.as_slice().into();
                record.push_aux(tag, Aux::ArrayU8(aux_array))?;
            }
            Self::ArrayI16(values) => {
                let aux_array: AuxArray<'_, i16> = values.as_slice().into();
                record.push_aux(tag, Aux::ArrayI16(aux_array))?;
            }
            Self::ArrayU16(values) => {
                let aux_array: AuxArray<'_, u16> = values.as_slice().into();
                record.push_aux(tag, Aux::ArrayU16(aux_array))?;
            }
            Self::ArrayI32(values) => {
                let aux_array: AuxArray<'_, i32> = values.as_slice().into();
                record.push_aux(tag, Aux::ArrayI32(aux_array))?;
            }
            Self::ArrayU32(values) => {
                let aux_array: AuxArray<'_, u32> = values.as_slice().into();
                record.push_aux(tag, Aux::ArrayU32(aux_array))?;
            }
            Self::ArrayFloat(bits) => {
                let values: Vec<f32> = bits.iter().copied().map(f32::from_bits).collect();
                let aux_array: AuxArray<'_, f32> = values.as_slice().into();
                record.push_aux(tag, Aux::ArrayFloat(aux_array))?;
            }
        }

        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct RepairResult {
    repaired: bool,
    report_name: bool,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct Summary {
    records_processed: u64,
    records_repaired: u64,
    report_names_written: u64,
}

#[derive(Debug)]
struct ProgressMeter {
    start: Instant,
    last_report: Instant,
    report_every: Option<Duration>,
    input_bytes: Option<u64>,
}

impl ProgressMeter {
    fn new(input_bytes: Option<u64>, progress_seconds: u64) -> Self {
        let now = Instant::now();
        Self {
            start: now,
            last_report: now,
            report_every: match progress_seconds {
                0 => None,
                seconds => Some(Duration::from_secs(seconds)),
            },
            input_bytes,
        }
    }

    fn maybe_report(&mut self, summary: &Summary, virtual_offset: i64) {
        let Some(report_every) = self.report_every else {
            return;
        };

        if self.last_report.elapsed() < report_every {
            return;
        }

        self.last_report = Instant::now();
        self.report("progress", summary, Some(virtual_offset));
    }

    fn finish(&self, summary: &Summary, virtual_offset: i64) {
        self.report("complete", summary, Some(virtual_offset));
    }

    fn report(&self, label: &str, summary: &Summary, virtual_offset: Option<i64>) {
        eprintln!(
            "[*] {label}: {}, {} repaired, {} reported, {}, {}, {} elapsed",
            format_count(summary.records_processed),
            format_count(summary.records_repaired),
            format_count(summary.report_names_written),
            self.format_percent_complete(virtual_offset),
            format_rate(summary.records_processed, self.start.elapsed()),
            format_duration(self.start.elapsed()),
        );
    }

    fn format_percent_complete(&self, virtual_offset: Option<i64>) -> String {
        let Some(input_bytes) = self.input_bytes.filter(|n| *n > 0) else {
            return "progress unavailable".to_string();
        };

        let Some(virtual_offset) = virtual_offset else {
            return "progress unavailable".to_string();
        };

        if virtual_offset < 0 {
            return "progress unavailable".to_string();
        }

        let compressed_offset = (virtual_offset as u64) >> 16;
        let percent = (compressed_offset as f64 * 100.0 / input_bytes as f64).min(100.0);
        format!("{percent:.1}% input")
    }
}

fn collect_target_values(record: &Record) -> Result<Vec<OwnedAux>> {
    let mut values = Vec::new();

    for result in record.aux_iter() {
        let (tag, value) = result.context("failed to parse auxiliary BAM field")?;
        if tag == TARGET_TAG {
            values.push(OwnedAux::from_aux(value));
        }
    }

    Ok(values)
}

fn repair_record(record: &mut Record) -> Result<RepairResult> {
    let mut mx_count = 0;
    for result in record.aux_iter() {
        let (tag, _) = result.context("failed to parse auxiliary BAM field")?;
        if tag == TARGET_TAG {
            mx_count += 1;
            if mx_count > 1 {
                break;
            }
        }
    }

    if mx_count <= 1 {
        return Ok(RepairResult::default());
    }

    let values = collect_target_values(record)?;

    let first = values
        .first()
        .cloned()
        .context("duplicate mx list unexpectedly had no first value")?;
    let has_mismatch = values.iter().skip(1).any(|value| value != &first);
    let report_name = values.len() > 2 || has_mismatch;

    while record.aux(TARGET_TAG).is_ok() {
        record.remove_aux(TARGET_TAG)?;
    }

    first.push_to_record(record, TARGET_TAG)?;

    Ok(RepairResult {
        repaired: true,
        report_name,
    })
}

fn read_name(record: &Record) -> Result<String> {
    let qname = record.qname();

    if qname.is_empty() {
        bail!("record needing mismatch reporting has an empty read name");
    }

    Ok(String::from_utf8_lossy(qname).into_owned())
}

fn run(
    input: &Path,
    output_bam: &Path,
    mismatches: &Path,
    threads: NonZeroUsize,
    progress_seconds: u64,
) -> Result<Summary> {
    let input_bytes = fs::metadata(input).ok().map(|metadata| metadata.len());
    let mut reader = bam::Reader::from_path(input)
        .with_context(|| format!("failed to open input BAM {}", input.display()))?;
    reader.set_threads(threads.get())?;
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(output_bam, &header, bam::Format::Bam)
        .with_context(|| format!("failed to create output BAM {}", output_bam.display()))?;
    writer.set_compression_level(CompressionLevel::Fastest)?;
    writer.set_threads(threads.get())?;
    let mismatch_file = File::create(mismatches)
        .with_context(|| format!("failed to create mismatch report {}", mismatches.display()))?;
    let mut mismatch_writer = BufWriter::new(mismatch_file);

    let mut summary = Summary::default();
    let mut progress = ProgressMeter::new(input_bytes, progress_seconds);

    let mut record = Record::new();

    while let Some(result) = reader.read(&mut record) {
        result.context("failed to read BAM record")?;
        summary.records_processed += 1;

        let repair = repair_record(&mut record).with_context(|| {
            format!(
                "failed to inspect/repair record {}",
                read_name_lossy(&record)
            )
        })?;

        if repair.repaired {
            summary.records_repaired += 1;
        }

        if repair.report_name {
            writeln!(mismatch_writer, "{}", read_name(&record)?)?;
            summary.report_names_written += 1;
        }

        writer
            .write(&record)
            .context("failed to write BAM record")?;

        progress.maybe_report(&summary, reader.tell());
    }

    mismatch_writer.flush()?;
    progress.finish(&summary, reader.tell());

    Ok(summary)
}

fn read_name_lossy(record: &Record) -> String {
    let qname = record.qname();

    if qname.is_empty() {
        "<empty-qname>".to_string()
    } else {
        String::from_utf8_lossy(qname).into_owned()
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let summary = run(
        &cli.input,
        &cli.output_bam,
        &cli.mismatches,
        cli.threads,
        cli.progress_seconds,
    )?;

    eprintln!(
        "[*] complete: {} records processed, {} records repaired, {} report names written",
        summary.records_processed, summary.records_repaired, summary.report_names_written
    );

    Ok(())
}

fn format_count(n: u64) -> String {
    match n {
        0..=999 => n.to_string(),
        1_000..=999_999 => format!("{:.1}k", n as f64 / 1_000.0),
        1_000_000..=999_999_999 => format!("{:.1}M", n as f64 / 1_000_000.0),
        _ => format!("{:.1}G", n as f64 / 1_000_000_000.0),
    }
}

fn format_duration(duration: Duration) -> String {
    let total_seconds = duration.as_secs();
    let hours = total_seconds / 3600;
    let minutes = (total_seconds % 3600) / 60;
    let seconds = total_seconds % 60;

    if hours > 0 {
        format!("{hours}:{minutes:02}:{seconds:02}")
    } else {
        format!("{minutes:02}:{seconds:02}")
    }
}

fn format_rate(records: u64, elapsed: Duration) -> String {
    let seconds = elapsed.as_secs_f64();

    if seconds <= 0.0 {
        return "0 rec/s".to_string();
    }

    let rate = (records as f64 / seconds).round() as u64;
    format!("{} rec/s", format_count(rate))
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::{fs, path::Path};

    use rust_htslib::bam::record::Record;
    use tempfile::tempdir;

    fn make_record(name: &[u8], mx_values: Vec<Aux<'_>>) -> Result<Record> {
        let tags: Vec<(&[u8], Aux<'_>)> = mx_values.into_iter().map(|v| (TARGET_TAG, v)).collect();
        make_record_complex(name, tags)
    }

    fn make_record_complex(name: &[u8], tags: Vec<(&[u8], Aux<'_>)>) -> Result<Record> {
        let mut record = Record::new();
        let qual = vec![30; 4];
        record.set(name, None, b"ACGT", &qual);

        for (tag, value) in tags {
            record.push_aux_unchecked(tag, value)?;
        }

        Ok(record)
    }

    fn write_input_bam(path: &Path, records: Vec<Record>) -> Result<()> {
        let header = bam::Header::new();
        let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam)?;

        for record in records {
            writer.write(&record)?;
        }

        Ok(())
    }

    fn read_records(path: &Path) -> Result<Vec<Record>> {
        let mut reader = bam::Reader::from_path(path)?;
        reader
            .records()
            .map(|result| result.map_err(anyhow::Error::from))
            .collect()
    }

    fn mx_values(record: &Record) -> Result<Vec<OwnedAux>> {
        collect_target_values(record)
    }

    fn all_tags(record: &Record) -> Result<Vec<(String, OwnedAux)>> {
        let mut tags = Vec::new();
        for result in record.aux_iter() {
            let (tag, value) = result.context("failed to parse auxiliary BAM field")?;
            let tag_str = String::from_utf8_lossy(tag).into_owned();
            tags.push((tag_str, OwnedAux::from_aux(value)));
        }
        Ok(tags)
    }

    #[test]
    fn repair_record_leaves_absent_or_single_mx_unchanged() -> Result<()> {
        let mut no_mx = make_record(b"no_mx", vec![])?;
        assert_eq!(repair_record(&mut no_mx)?, RepairResult::default());
        assert!(mx_values(&no_mx)?.is_empty());

        let mut one_mx = make_record(b"one_mx", vec![Aux::String("pass")])?;
        assert_eq!(repair_record(&mut one_mx)?, RepairResult::default());
        assert_eq!(
            mx_values(&one_mx)?,
            vec![OwnedAux::String("pass".to_string())]
        );

        Ok(())
    }

    #[test]
    fn repair_record_drops_identical_duplicate_without_report() -> Result<()> {
        let mut record = make_record(b"same", vec![Aux::String("pass"), Aux::String("pass")])?;

        assert_eq!(
            repair_record(&mut record)?,
            RepairResult {
                repaired: true,
                report_name: false
            }
        );
        assert_eq!(
            mx_values(&record)?,
            vec![OwnedAux::String("pass".to_string())]
        );

        Ok(())
    }

    #[test]
    fn repair_record_drops_mismatched_duplicate_and_reports() -> Result<()> {
        let mut record = make_record(b"different", vec![Aux::String("pass"), Aux::String("fail")])?;

        assert_eq!(
            repair_record(&mut record)?,
            RepairResult {
                repaired: true,
                report_name: true
            }
        );
        assert_eq!(
            mx_values(&record)?,
            vec![OwnedAux::String("pass".to_string())]
        );

        Ok(())
    }

    #[test]
    fn repair_record_reports_three_or_more_even_when_equal() -> Result<()> {
        let mut record = make_record(b"three", vec![Aux::I32(7), Aux::I32(7), Aux::I32(7)])?;

        assert_eq!(
            repair_record(&mut record)?,
            RepairResult {
                repaired: true,
                report_name: true
            }
        );
        assert_eq!(mx_values(&record)?, vec![OwnedAux::I32(7)]);

        Ok(())
    }

    #[test]
    fn repair_record_preserves_float_bits_for_exact_comparison() -> Result<()> {
        let nan = f32::from_bits(0x7fc0_0001);
        let mut record = make_record(b"nan", vec![Aux::Float(nan), Aux::Float(nan)])?;

        assert_eq!(
            repair_record(&mut record)?,
            RepairResult {
                repaired: true,
                report_name: false
            }
        );
        assert_eq!(mx_values(&record)?, vec![OwnedAux::Float(nan.to_bits())]);

        Ok(())
    }

    #[test]
    fn run_repairs_bam_and_writes_mismatch_names() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("input.bam");
        let output = dir.path().join("output.bam");
        let mismatches = dir.path().join("mismatches.txt");

        write_input_bam(
            &input,
            vec![
                make_record(b"same", vec![Aux::String("x"), Aux::String("x")])?,
                make_record(b"different", vec![Aux::String("x"), Aux::String("y")])?,
                make_record(b"three", vec![Aux::U8(1), Aux::U8(1), Aux::U8(1)])?,
                make_record(b"single", vec![Aux::String("z")])?,
            ],
        )?;

        let summary = run(&input, &output, &mismatches, NonZeroUsize::MIN, 0)?;

        assert_eq!(
            summary,
            Summary {
                records_processed: 4,
                records_repaired: 3,
                report_names_written: 2
            }
        );

        let records = read_records(&output)?;
        assert_eq!(records.len(), 4);
        assert_eq!(
            records.iter().map(mx_values).collect::<Result<Vec<_>>>()?,
            vec![
                vec![OwnedAux::String("x".to_string())],
                vec![OwnedAux::String("x".to_string())],
                vec![OwnedAux::U8(1)],
                vec![OwnedAux::String("z".to_string())],
            ]
        );

        assert_eq!(fs::read_to_string(mismatches)?, "different\nthree\n");

        Ok(())
    }

    #[test]
    fn repair_record_preserves_surrounding_tags() -> Result<()> {
        let mut record = make_record_complex(
            b"interleaved",
            vec![
                (b"NM", Aux::I32(1)),
                (TARGET_TAG, Aux::String("v1")),
                (b"AS", Aux::I32(30)),
                (TARGET_TAG, Aux::String("v1")),
            ],
        )?;

        repair_record(&mut record)?;

        let tags = all_tags(&record)?;
        // Order changes because of remove/push logic, but tags must be present
        assert!(tags.contains(&("NM".to_string(), OwnedAux::I32(1))));
        assert!(tags.contains(&("AS".to_string(), OwnedAux::I32(30))));
        assert!(tags.contains(&("mx".to_string(), OwnedAux::String("v1".to_string()))));
        assert_eq!(tags.len(), 3);

        Ok(())
    }

    #[test]
    fn repair_record_reports_heterogeneous_type_mismatch() -> Result<()> {
        let mut record = make_record(b"types", vec![Aux::I32(10), Aux::Float(10.0)])?;

        let res = repair_record(&mut record)?;
        assert!(res.report_name);
        assert_eq!(mx_values(&record)?, vec![OwnedAux::I32(10)]);

        Ok(())
    }

    #[test]
    fn repair_record_handles_empty_collections() -> Result<()> {
        let mut record = make_record(b"empty", vec![Aux::String(""), Aux::String("")])?;
        let res = repair_record(&mut record)?;
        assert!(!res.report_name);
        assert_eq!(mx_values(&record)?, vec![OwnedAux::String("".to_string())]);

        Ok(())
    }

    #[test]
    fn read_name_handles_non_utf8_gracefully() -> Result<()> {
        let invalid_utf8 = b"\xff\xfe\xfd";
        let record = make_record(invalid_utf8, vec![])?;

        // Lossy conversion replaces invalid bytes with replacement characters (U+FFFD)
        let expected = "\u{FFFD}\u{FFFD}\u{FFFD}";
        assert_eq!(read_name_lossy(&record), expected);
        assert_eq!(read_name(&record)?, expected);

        Ok(())
    }

    #[test]
    fn read_name_handles_max_length_qname() -> Result<()> {
        let max_name = vec![b'A'; 254];
        let record = make_record(&max_name, vec![])?;

        assert_eq!(read_name(&record)?.len(), 254);

        Ok(())
    }

    #[test]
    fn cli_rejects_zero_threads() {
        let err = Cli::try_parse_from([
            "lr_drop_duptag",
            "--input",
            "in.bam",
            "--output-bam",
            "out.bam",
            "--mismatches",
            "mismatches.txt",
            "--threads",
            "0",
        ])
        .expect_err("zero threads must be rejected");

        assert!(err.to_string().contains("invalid value"));
    }

    #[test]
    fn progress_seconds_zero_disables_periodic_reports() {
        let meter = ProgressMeter::new(Some(100), 0);
        assert!(meter.report_every.is_none());
    }

    #[test]
    fn progress_meter_formats_bgzf_virtual_offset_percent() {
        let meter = ProgressMeter::new(Some(1_000), 60);
        let virtual_offset = 250_i64 << 16;

        assert_eq!(
            meter.format_percent_complete(Some(virtual_offset)),
            "25.0% input"
        );
    }

    #[test]
    fn progress_meter_handles_unknown_progress() {
        let meter = ProgressMeter::new(None, 60);

        assert_eq!(
            meter.format_percent_complete(Some(250_i64 << 16)),
            "progress unavailable"
        );
        assert_eq!(
            meter.format_percent_complete(Some(-1)),
            "progress unavailable"
        );
    }
}
