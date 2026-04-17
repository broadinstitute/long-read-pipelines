use std::env;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::sync::{
    Arc,
    atomic::{AtomicU64, Ordering},
};
use std::time::{Duration, Instant};

use bio::alphabets::dna;
use flate2::{Compression, read::MultiGzDecoder, write::GzEncoder};
use rustc_hash::FxHashMap;
use noodles::bam;
use noodles::fastq;
use noodles::sam::{
    alignment::{
        io::Write as _,
        record::{Cigar as _, QualityScores as _},
        record::cigar::{op::Kind, Op},
        record_buf::{Cigar, QualityScores, Sequence},
        RecordBuf,
    },
};

type FastqCache = FxHashMap<Vec<u8>, (Vec<u8>, Vec<u8>)>;
const REPORT_EVERY: Duration = Duration::from_secs(60);

#[derive(Clone, Debug, Default)]
struct ByteCounter(Arc<AtomicU64>);

impl ByteCounter {
    fn bytes_read(&self) -> u64 {
        self.0.load(Ordering::Relaxed)
    }

    fn add(&self, n: u64) {
        self.0.fetch_add(n, Ordering::Relaxed);
    }
}

#[derive(Debug)]
struct CountingReader<R> {
    inner: R,
    counter: ByteCounter,
}

impl<R> CountingReader<R> {
    fn new(inner: R) -> (Self, ByteCounter) {
        let counter = ByteCounter::default();
        (
            Self {
                inner,
                counter: counter.clone(),
            },
            counter,
        )
    }
}

impl<R> Read for CountingReader<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let bytes_read = self.inner.read(buf)?;
        self.counter.add(bytes_read as u64);
        Ok(bytes_read)
    }
}

#[derive(Debug)]
struct ProgressMeter {
    phase: &'static str,
    start: Instant,
    last_report: Instant,
    report_every: Duration,
    total_bytes: Option<u64>,
}

impl ProgressMeter {
    fn new(phase: &'static str, total_bytes: Option<u64>) -> Self {
        let now = Instant::now();
        Self {
            phase,
            start: now,
            last_report: now,
            report_every: REPORT_EVERY,
            total_bytes,
        }
    }

    fn maybe_report_fastq(&mut self, bytes_read: u64, reads: u64) {
        if !self.should_report() {
            return;
        }

        self.last_report = Instant::now();
        eprintln!(
            "[*] {}: {}, {}, {}, {}, {} elapsed",
            self.phase,
            format_count(reads),
            self.format_percent_complete(bytes_read),
            format_rate(reads, self.start.elapsed(), "reads/s"),
            self.format_rss(),
            format_duration(self.start.elapsed()),
        );
    }

    fn maybe_report_bam(&mut self, bytes_read: Option<u64>, records: u64, patched: u64) {
        if !self.should_report() {
            return;
        }

        self.last_report = Instant::now();
        eprintln!(
            "[*] {}: {}, {} patched, {}, {}, {}, {}, {} elapsed",
            self.phase,
            format_count(records),
            format_count(patched),
            format_fraction(patched, records, "patched"),
            self.format_optional_percent_complete(bytes_read),
            format_rate(records, self.start.elapsed(), "rec/s"),
            self.format_rss(),
            format_duration(self.start.elapsed()),
        );
    }

    fn finish_fastq(&self, bytes_read: u64, reads: u64) {
        eprintln!(
            "[*] {} complete: {}, {}, {}, {}, {} elapsed",
            self.phase,
            format_count(reads),
            self.format_percent_complete(bytes_read),
            format_rate(reads, self.start.elapsed(), "reads/s"),
            self.format_rss(),
            format_duration(self.start.elapsed()),
        );
    }

    fn finish_bam(&self, bytes_read: Option<u64>, records: u64, patched: u64) {
        eprintln!(
            "[*] {} complete: {}, {} patched, {}, {}, {}, {}, {} elapsed",
            self.phase,
            format_count(records),
            format_count(patched),
            format_fraction(patched, records, "patched"),
            self.format_optional_percent_complete(bytes_read),
            format_rate(records, self.start.elapsed(), "rec/s"),
            self.format_rss(),
            format_duration(self.start.elapsed()),
        );
    }

    fn should_report(&self) -> bool {
        self.last_report.elapsed() >= self.report_every
    }

    fn format_percent_complete(&self, bytes_read: u64) -> String {
        match self.total_bytes {
            Some(total_bytes) if total_bytes > 0 => {
                let percent = (bytes_read as f64 * 100.0 / total_bytes as f64).min(100.0);
                format!("{percent:.1}% complete")
            }
            _ => "progress unavailable".to_string(),
        }
    }

    fn format_optional_percent_complete(&self, bytes_read: Option<u64>) -> String {
        match bytes_read {
            Some(bytes_read) => self.format_percent_complete(bytes_read),
            None => "progress unavailable".to_string(),
        }
    }

    fn format_rss(&self) -> String {
        current_rss_bytes()
            .map(|rss_bytes| format!("RSS {}", format_bytes_gib(rss_bytes)))
            .unwrap_or_else(|| "RSS unavailable".to_string())
    }
}

fn format_bytes_gib(bytes: u64) -> String {
    format!("{:.1} GiB", bytes as f64 / 1024.0 / 1024.0 / 1024.0)
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

fn format_rate(units: u64, elapsed: Duration, suffix: &str) -> String {
    let secs = elapsed.as_secs_f64();
    if secs <= 0.0 {
        return format!("0 {suffix}");
    }

    let rate = (units as f64 / secs).round() as u64;
    format!("{} {}", format_count(rate), suffix)
}

fn format_fraction(numerator: u64, denominator: u64, suffix: &str) -> String {
    if denominator == 0 {
        format!("0.0% {suffix}")
    } else {
        format!("{:.1}% {}", numerator as f64 * 100.0 / denominator as f64, suffix)
    }
}

fn current_rss_bytes() -> Option<u64> {
    let status = fs::read_to_string("/proc/self/status").ok()?;

    status
        .lines()
        .find(|line| line.starts_with("VmRSS:"))
        .and_then(|line| line.split_whitespace().nth(1))
        .and_then(|kb| kb.parse::<u64>().ok())
        .map(|kb| kb * 1024)
}

fn open_fastq_reader(
    fastq_path: &str,
) -> anyhow::Result<(fastq::Reader<BufReader<Box<dyn Read>>>, ByteCounter, Option<u64>)> {
    let fastq_total_bytes = fs::metadata(fastq_path).ok().map(|m| m.len());
    let fastq_file = File::open(fastq_path)?;
    let (counting_fastq_reader, fastq_counter) = CountingReader::new(fastq_file);
    let mut raw_fastq_reader = BufReader::new(counting_fastq_reader);
    let is_gzip = raw_fastq_reader.fill_buf()?.starts_with(&[0x1f, 0x8b]);

    let decoded_reader: Box<dyn Read> = if is_gzip {
        Box::new(MultiGzDecoder::new(raw_fastq_reader))
    } else {
        Box::new(raw_fastq_reader)
    };

    Ok((
        fastq::Reader::new(BufReader::new(decoded_reader)),
        fastq_counter,
        fastq_total_bytes,
    ))
}

fn read_consuming_cigar_len(cigar: &Cigar) -> anyhow::Result<usize> {
    cigar.iter().try_fold(0usize, |len, op_res| {
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

fn record_type(record: &RecordBuf) -> &'static str {
    let flags = record.flags();

    if flags.is_secondary() {
        "secondary"
    } else if flags.is_supplementary() {
        "supplementary"
    } else {
        "primary"
    }
}

fn has_hard_clip(cigar: &Cigar) -> anyhow::Result<bool> {
    for op_res in cigar.iter() {
        let op = op_res?;
        if op.kind() == Kind::HardClip {
            return Ok(true);
        }
    }

    Ok(false)
}

fn replace_hard_clips_with_soft_clips(cigar: &Cigar) -> anyhow::Result<Cigar> {
    let mut new_ops: Vec<Op> = Vec::new();

    for op_res in cigar.iter() {
        let op = op_res?;
        let kind = if op.kind() == Kind::HardClip {
            Kind::SoftClip
        } else {
            op.kind()
        };

        if let Some(last_op) = new_ops.last_mut() {
            if last_op.kind() == kind {
                *last_op = Op::new(kind, last_op.len() + op.len());
                continue;
            }
        }

        new_ops.push(Op::new(kind, op.len()));
    }

    Ok(Cigar::from(new_ops))
}

fn decode_fastq_quality_scores(record: &fastq::Record) -> anyhow::Result<Vec<u8>> {
    record
        .quality_scores()
        .iter()
        .map(|&q| {
            anyhow::ensure!(
                q >= 33,
                "invalid FASTQ quality byte {} for read {}",
                q,
                String::from_utf8_lossy(record.name().as_ref()),
            );
            Ok(q - 33)
        })
        .collect()
}

#[cfg(test)]
fn load_fastq_cache<R>(reader: R) -> anyhow::Result<(FastqCache, usize)>
where
    R: BufRead,
{
    let mut fastq_reader = fastq::Reader::new(reader);
    let mut cache = FxHashMap::default();
    let mut read_count = 0;

    for result in fastq_reader.records() {
        let record = result?;
        let quality_scores = decode_fastq_quality_scores(&record)?;
        // We store Name -> (Sequence, Qualities)
        cache.insert(
            record.name().to_vec(),
            (record.sequence().to_vec(), quality_scores),
        );
        read_count += 1;
        if read_count % 1_000_000 == 0 {
            println!("    Loaded {}M reads...", read_count / 1_000_000);
        }
    }

    Ok((cache, read_count))
}

fn restore_record(record: &mut RecordBuf, cache: &FastqCache) -> anyhow::Result<bool> {
    let cigar = record.cigar().clone();

    if !has_hard_clip(&cigar)? {
        return Ok(false);
    }

    let read_name_bytes = record
        .name()
        .ok_or_else(|| {
            anyhow::anyhow!(
                "hard-clipped {} alignment is missing a read name",
                record_type(record)
            )
        })?
        .to_vec();

    let (raw_seq, raw_qual) = cache.get(read_name_bytes.as_slice()).ok_or_else(|| {
        anyhow::anyhow!(
            "missing FASTQ entry for hard-clipped {} alignment {}",
            record_type(record),
            String::from_utf8_lossy(&read_name_bytes),
        )
    })?;

    let mut final_seq = raw_seq.clone();
    let mut final_qual = raw_qual.clone();
    let fastq_seq_len = final_seq.len();
    let fastq_qual_len = final_qual.len();

    if record.flags().is_reverse_complemented() {
        final_seq = dna::revcomp(&final_seq);
        final_qual.reverse();
    }

    *record.sequence_mut() = Sequence::from(final_seq);
    *record.quality_scores_mut() = QualityScores::from(final_qual);
    *record.cigar_mut() = replace_hard_clips_with_soft_clips(&cigar)?;

    let restored_seq_len = record.sequence().len();
    let restored_qual_len = record.quality_scores().len();
    let cigar_len = read_consuming_cigar_len(record.cigar())?;

    anyhow::ensure!(
        fastq_seq_len == fastq_qual_len
            && fastq_seq_len == restored_seq_len
            && fastq_seq_len == restored_qual_len
            && fastq_seq_len == cigar_len,
        "length mismatch for {} alignment {}: FASTQ seq={}, FASTQ qual={}, restored seq={}, restored qual={}, CIGAR read span={}",
        record_type(record),
        String::from_utf8_lossy(&read_name_bytes),
        fastq_seq_len,
        fastq_qual_len,
        restored_seq_len,
        restored_qual_len,
        cigar_len,
    );

    Ok(true)
}

fn run(bam_path: &str, fastq_path: &str, out_path: &str) -> anyhow::Result<usize> {
    // --- PASS 1: Load FASTQ into RAM ---
    // With 96GB RAM, we cache everything using FxHashMap for O(1) lookups.
    eprintln!("[*] Loading FASTQ into RAM cache...");
    let (mut fastq_reader, fastq_counter, fastq_total_bytes) = open_fastq_reader(fastq_path)?;
    let mut cache = FxHashMap::default();
    let mut read_count = 0u64;
    let mut fastq_meter = ProgressMeter::new("FASTQ cache", fastq_total_bytes);

    for result in fastq_reader.records() {
        let record = result?;
        let quality_scores = decode_fastq_quality_scores(&record)?;
        cache.insert(
            record.name().to_vec(),
            (record.sequence().to_vec(), quality_scores),
        );
        read_count += 1;
        fastq_meter.maybe_report_fastq(fastq_counter.bytes_read(), read_count);
    }
    fastq_meter.finish_fastq(fastq_counter.bytes_read(), read_count);

    // --- PASS 2: Process BAM and Patch Hard Clips ---
    eprintln!("[*] Opening BAM for patching...");
    let bam_total_bytes = fs::metadata(bam_path).ok().map(|m| m.len());
    let bam_file = File::open(bam_path)?;
    let (counting_bam_reader, bam_counter) = CountingReader::new(bam_file);
    let mut reader = bam::io::reader::Builder::default().build_from_reader(counting_bam_reader);
    let header = reader.read_header()?;
    
    let mut writer = bam::io::writer::Builder::default().build_from_path(out_path)?;
    writer.write_header(&header)?;

    let mut patched_count = 0u64;
    let mut records_processed = 0u64;
    let mut bam_meter = ProgressMeter::new("BAM patching", bam_total_bytes);

    for result in reader.record_bufs(&header) {
        let mut record = result?;
        records_processed += 1;
        if restore_record(&mut record, &cache)? {
            patched_count += 1;
        }
        writer.write_alignment_record(&header, &record)?;
        bam_meter.maybe_report_bam(Some(bam_counter.bytes_read()), records_processed, patched_count);
    }

    writer.try_finish()?;
    bam_meter.finish_bam(Some(bam_counter.bytes_read()), records_processed, patched_count);

    Ok(patched_count as usize)
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!("Usage: bam_restorer <in.bam> <in.fastq> <out.bam>");
        std::process::exit(1);
    }

    let patched_count = run(&args[1], &args[2], &args[3])?;
    eprintln!("[*] Success! Patched {} alignments.", patched_count);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::fs;
    use std::io::{Cursor, Write};
    use std::num::NonZeroUsize;
    use std::path::{Path, PathBuf};
    use std::time::{SystemTime, UNIX_EPOCH};

    use noodles::core::Position;
    use noodles::sam::{
        self,
        alignment::record::{Flags, MappingQuality},
    };

    fn cigar(ops: &[(Kind, usize)]) -> Cigar {
        ops.iter()
            .map(|(kind, len)| Op::new(*kind, *len))
            .collect()
    }

    fn make_cache(entries: &[(&[u8], &[u8], &[u8])]) -> FastqCache {
        let mut cache = FastqCache::default();

        for (name, seq, qual) in entries {
            cache.insert(name.to_vec(), (seq.to_vec(), qual.to_vec()));
        }

        cache
    }

    fn make_record(name: Option<&str>, flags: Flags, cigar: Cigar, seq: &[u8], qual: &[u8]) -> RecordBuf {
        let mut builder = RecordBuf::builder()
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::MIN)
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_cigar(cigar)
            .set_sequence(Sequence::from(seq))
            .set_quality_scores(QualityScores::from(qual.to_vec()));

        if let Some(name) = name {
            builder = builder.set_name(name);
        }

        builder.build()
    }

    fn temp_dir(prefix: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "bam_restorer_{prefix}_{}_{}",
            std::process::id(),
            unique
        ));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write_test_fastq(path: &Path, records: &[(&str, &str, &str)]) {
        let mut file = fs::File::create(path).unwrap();

        for (name, seq, qual) in records {
            writeln!(file, "@{name}").unwrap();
            writeln!(file, "{seq}").unwrap();
            writeln!(file, "+").unwrap();
            writeln!(file, "{qual}").unwrap();
        }
    }

    fn write_test_fastq_gz(path: &Path, records: &[(&str, &str, &str)]) {
        let file = fs::File::create(path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());

        for (name, seq, qual) in records {
            writeln!(encoder, "@{name}").unwrap();
            writeln!(encoder, "{seq}").unwrap();
            writeln!(encoder, "+").unwrap();
            writeln!(encoder, "{qual}").unwrap();
        }

        encoder.finish().unwrap();
    }

    fn write_test_bam(path: &Path, records: &[RecordBuf]) -> sam::Header {
        let header = sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                sam::header::record::value::Map::<sam::header::record::value::map::ReferenceSequence>::new(
                    NonZeroUsize::new(100_000).unwrap(),
                ),
            )
            .build();

        let mut writer = bam::io::writer::Builder::default()
            .build_from_path(path)
            .unwrap();
        writer.write_header(&header).unwrap();

        for record in records {
            writer.write_alignment_record(&header, record).unwrap();
        }

        writer.try_finish().unwrap();

        header
    }

    fn read_bam_records(path: &Path) -> (sam::Header, Vec<RecordBuf>) {
        let mut reader = bam::io::reader::Builder::default()
            .build_from_path(path)
            .unwrap();
        let header = reader.read_header().unwrap();
        let records = reader
            .record_bufs(&header)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        (header, records)
    }

    fn write_fastq_from_primary_records(input_bam_path: &Path, fastq_path: &Path) {
        let (_header, records) = read_bam_records(input_bam_path);
        let mut file = fs::File::create(fastq_path).unwrap();

        for record in records {
            if record.flags().is_secondary() || record.flags().is_supplementary() {
                continue;
            }

            let mut seq = record.sequence().as_ref().to_vec();
            let mut qual = record.quality_scores().as_ref().to_vec();

            if record.flags().is_reverse_complemented() {
                seq = dna::revcomp(&seq);
                qual.reverse();
            }

            let qual_ascii: Vec<u8> = qual.into_iter().map(|q| q + 33).collect();

            file.write_all(b"@").unwrap();
            file.write_all(record.name().unwrap().as_ref()).unwrap();
            file.write_all(b"\n").unwrap();
            file.write_all(&seq).unwrap();
            file.write_all(b"\n+\n").unwrap();
            file.write_all(&qual_ascii).unwrap();
            file.write_all(b"\n").unwrap();
        }
    }

    #[test]
    fn read_consuming_cigar_len_counts_only_read_consuming_ops() {
        let cigar = cigar(&[
            (Kind::HardClip, 5),
            (Kind::SoftClip, 10),
            (Kind::Match, 20),
            (Kind::Insertion, 2),
            (Kind::Deletion, 3),
            (Kind::Skip, 4),
            (Kind::Pad, 5),
            (Kind::SequenceMatch, 6),
            (Kind::SequenceMismatch, 7),
            (Kind::HardClip, 8),
        ]);

        assert_eq!(read_consuming_cigar_len(&cigar).unwrap(), 45);
    }

    #[test]
    fn has_hard_clip_detects_presence_and_absence() {
        assert!(has_hard_clip(&cigar(&[(Kind::HardClip, 3), (Kind::Match, 4)])).unwrap());
        assert!(!has_hard_clip(&cigar(&[(Kind::SoftClip, 3), (Kind::Match, 4)])).unwrap());
    }

    #[test]
    fn replace_hard_clips_with_soft_clips_preserves_other_ops() {
        let original = cigar(&[
            (Kind::HardClip, 2),
            (Kind::Match, 4),
            (Kind::Insertion, 1),
            (Kind::HardClip, 3),
        ]);
        let replaced = replace_hard_clips_with_soft_clips(&original).unwrap();
        let ops = replaced
            .iter()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(ops[0], Op::new(Kind::SoftClip, 2));
        assert_eq!(ops[1], Op::new(Kind::Match, 4));
        assert_eq!(ops[2], Op::new(Kind::Insertion, 1));
        assert_eq!(ops[3], Op::new(Kind::SoftClip, 3));
    }

    #[test]
    fn replace_hard_clips_with_soft_clips_merges_adjacent_soft_clips() {
        let original = cigar(&[
            (Kind::SoftClip, 10),
            (Kind::HardClip, 20),
            (Kind::Match, 100),
            (Kind::HardClip, 30),
            (Kind::SoftClip, 5),
        ]);
        let replaced = replace_hard_clips_with_soft_clips(&original).unwrap();
        let ops = replaced
            .iter()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0], Op::new(Kind::SoftClip, 30));
        assert_eq!(ops[1], Op::new(Kind::Match, 100));
        assert_eq!(ops[2], Op::new(Kind::SoftClip, 35));
    }

    #[test]
    fn record_type_classifies_primary_secondary_and_supplementary() {
        assert_eq!(record_type(&make_record(Some("r"), Flags::empty(), cigar(&[(Kind::Match, 4)]), b"ACGT", b"IIII")), "primary");
        assert_eq!(record_type(&make_record(Some("r"), Flags::SECONDARY, cigar(&[(Kind::Match, 4)]), b"ACGT", b"IIII")), "secondary");
        assert_eq!(record_type(&make_record(Some("r"), Flags::SUPPLEMENTARY, cigar(&[(Kind::Match, 4)]), b"ACGT", b"IIII")), "supplementary");
    }

    #[test]
    fn load_fastq_cache_counts_records_and_last_duplicate_wins() {
        let fastq = b"@r1\nACGT\n+\n!!!!\n@r2\nTGCA\n+\n####\n@r1\nNNNN\n+\n$$$$\n";
        let (cache, read_count) = load_fastq_cache(Cursor::new(&fastq[..])).unwrap();

        assert_eq!(read_count, 3);
        assert_eq!(cache.len(), 2);
        assert_eq!(cache.get(b"r1".as_slice()).unwrap().0, b"NNNN");
        assert_eq!(cache.get(b"r1".as_slice()).unwrap().1, vec![3, 3, 3, 3]);
    }

    #[test]
    fn decode_fastq_quality_scores_rejects_invalid_ascii() {
        let record = fastq::Record::new(
            fastq::record::Definition::new("bad", ""),
            "ACGT",
            vec![32, 33, 34, 35],
        );
        let err = decode_fastq_quality_scores(&record).unwrap_err();

        assert!(err.to_string().contains("invalid FASTQ quality byte 32"));
    }

    #[test]
    fn restore_record_ignores_unclipped_record_without_name() {
        let original = make_record(None, Flags::empty(), cigar(&[(Kind::Match, 4)]), b"ACGT", b"IIII");
        let mut record = original.clone();
        let cache = FastqCache::default();

        assert!(!restore_record(&mut record, &cache).unwrap());
        assert_eq!(record, original);
    }

    #[test]
    fn restore_record_restores_forward_hardclipped_record() {
        let mut record = make_record(
            Some("read1"),
            Flags::SUPPLEMENTARY,
            cigar(&[(Kind::HardClip, 2), (Kind::Match, 5)]),
            b"GTRYN",
            b"JKLMN",
        );
        let cache = make_cache(&[(b"read1", b"ACGTRYN", &[32, 33, 34, 35, 36, 37, 38])]);

        assert!(restore_record(&mut record, &cache).unwrap());
        assert_eq!(record.sequence().as_ref(), b"ACGTRYN");
        assert_eq!(record.quality_scores().as_ref(), &[32, 33, 34, 35, 36, 37, 38]);
        assert_eq!(read_consuming_cigar_len(record.cigar()).unwrap(), 7);

        let ops = record.cigar().iter().collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(ops, vec![Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 5)]);
    }

    #[test]
    fn restore_record_restores_reverse_record_and_reverses_qualities_with_iupac() {
        let mut record = make_record(
            Some("read2"),
            Flags::SECONDARY | Flags::REVERSE_COMPLEMENTED,
            cigar(&[(Kind::HardClip, 2), (Kind::Match, 5)]),
            b"NNNNN",
            b"!!!!!",
        );
        let cache = make_cache(&[(b"read2", b"ACGTRYN", &[32, 33, 34, 35, 36, 37, 38])]);

        assert!(restore_record(&mut record, &cache).unwrap());
        assert_eq!(record.sequence().as_ref(), dna::revcomp(b"ACGTRYN"));
        assert_eq!(record.quality_scores().as_ref(), &[38, 37, 36, 35, 34, 33, 32]);
    }

    #[test]
    fn restore_record_restores_reverse_record_with_lowercase_and_iupac() {
        let mut record = make_record(
            Some("read2"),
            Flags::SECONDARY | Flags::REVERSE_COMPLEMENTED,
            cigar(&[(Kind::HardClip, 2), (Kind::Match, 5)]),
            b"NNNNN",
            b"!!!!!",
        );
        let cache = make_cache(&[(b"read2", b"acgTRyn", &[32, 33, 34, 35, 36, 37, 38])]);

        assert!(restore_record(&mut record, &cache).unwrap());
        assert_eq!(record.sequence().as_ref(), b"nrYAcgt");
        assert_eq!(record.quality_scores().as_ref(), &[38, 37, 36, 35, 34, 33, 32]);
    }

    #[test]
    fn restore_record_errors_for_hardclipped_record_without_name() {
        let mut record = make_record(None, Flags::SUPPLEMENTARY, cigar(&[(Kind::HardClip, 2), (Kind::Match, 5)]), b"AAAAA", b"!!!!!");
        let err = restore_record(&mut record, &FastqCache::default()).unwrap_err();

        assert!(err.to_string().contains("missing a read name"));
        assert!(err.to_string().contains("supplementary alignment"));
    }

    #[test]
    fn restore_record_errors_for_missing_cache_entry() {
        let mut record = make_record(
            Some("read3"),
            Flags::SECONDARY,
            cigar(&[(Kind::HardClip, 1), (Kind::Match, 6)]),
            b"CCCCCC",
            b"!!!!!!",
        );
        let err = restore_record(&mut record, &FastqCache::default()).unwrap_err();

        assert!(err.to_string().contains("missing FASTQ entry"));
        assert!(err.to_string().contains("secondary alignment read3"));
    }

    #[test]
    fn restore_record_errors_for_fastq_sequence_quality_length_mismatch() {
        let mut record = make_record(
            Some("read4"),
            Flags::SUPPLEMENTARY,
            cigar(&[(Kind::HardClip, 2), (Kind::Match, 5)]),
            b"TTTTT",
            b"!!!!!",
        );
        let cache = make_cache(&[(b"read4", b"ACGTRYN", &[32, 33, 34, 35, 36, 37])]);
        let err = restore_record(&mut record, &cache).unwrap_err();

        assert!(err.to_string().contains("length mismatch"));
        assert!(err.to_string().contains("FASTQ seq=7, FASTQ qual=6"));
    }

    #[test]
    fn restore_record_errors_for_cigar_length_mismatch() {
        let mut record = make_record(
            Some("read5"),
            Flags::SUPPLEMENTARY,
            cigar(&[(Kind::HardClip, 2), (Kind::Match, 3)]),
            b"AAA",
            b"!!!",
        );
        let cache = make_cache(&[(b"read5", b"ACGTRYN", &[32, 33, 34, 35, 36, 37, 38])]);
        let err = restore_record(&mut record, &cache).unwrap_err();

        assert!(err.to_string().contains("length mismatch"));
        assert!(err.to_string().contains("CIGAR read span=5"));
    }

    #[test]
    fn run_restores_hardclipped_records_end_to_end() {
        let dir = temp_dir("run");
        let fastq_path = dir.join("reads.fastq");
        let input_bam_path = dir.join("input.bam");
        let output_bam_path = dir.join("output.bam");

        write_test_fastq(
            &fastq_path,
            &[
                ("read1", "ACGTRYN", "ABCDEFG"),
                ("read2", "TTACCGG", "HIJKLMN"),
            ],
        );

        let primary = make_record(
            Some("read1"),
            Flags::empty(),
            cigar(&[(Kind::Match, 7)]),
            b"ACGTRYN",
            &[32, 33, 34, 35, 36, 37, 38],
        );
        let secondary = make_record(
            Some("read1"),
            Flags::SECONDARY,
            cigar(&[(Kind::HardClip, 1), (Kind::Match, 6)]),
            b"CGTRYN",
            &[33, 34, 35, 36, 37, 38],
        );
        let supplementary = make_record(
            Some("read2"),
            Flags::SUPPLEMENTARY | Flags::REVERSE_COMPLEMENTED,
            cigar(&[(Kind::HardClip, 2), (Kind::Match, 5)]),
            b"CCGGT",
            &[39, 38, 37, 36, 35],
        );

        write_test_bam(&input_bam_path, &[primary, secondary, supplementary]);

        let patched_count = run(
            input_bam_path.to_str().unwrap(),
            fastq_path.to_str().unwrap(),
            output_bam_path.to_str().unwrap(),
        )
        .unwrap();

        assert_eq!(patched_count, 2);

        let (_header, records) = read_bam_records(&output_bam_path);
        assert_eq!(records.len(), 3);

        let read1_primary = records
            .iter()
            .find(|r| {
                let name: &[u8] = r.name().unwrap().as_ref();
                name == b"read1" && !r.flags().is_secondary()
            })
            .unwrap();
        assert_eq!(read1_primary.sequence().as_ref(), b"ACGTRYN");
        assert_eq!(read1_primary.quality_scores().as_ref(), &[32, 33, 34, 35, 36, 37, 38]);
        assert_eq!(
            read1_primary
                .cigar()
                .iter()
                .collect::<Result<Vec<_>, _>>()
                .unwrap(),
            vec![Op::new(Kind::Match, 7)]
        );

        let read1_secondary = records
            .iter()
            .find(|r| {
                let name: &[u8] = r.name().unwrap().as_ref();
                name == b"read1" && r.flags().is_secondary()
            })
            .unwrap();
        assert_eq!(read1_secondary.sequence().as_ref(), b"ACGTRYN");
        assert_eq!(read1_secondary.quality_scores().as_ref(), &[32, 33, 34, 35, 36, 37, 38]);
        assert_eq!(
            read1_secondary
                .cigar()
                .iter()
                .collect::<Result<Vec<_>, _>>()
                .unwrap(),
            vec![Op::new(Kind::SoftClip, 1), Op::new(Kind::Match, 6)]
        );

        let read2_supplementary = records
            .iter()
            .find(|r| {
                let name: &[u8] = r.name().unwrap().as_ref();
                name == b"read2"
            })
            .unwrap();
        assert_eq!(read2_supplementary.sequence().as_ref(), dna::revcomp(b"TTACCGG"));
        assert_eq!(read2_supplementary.quality_scores().as_ref(), &[45, 44, 43, 42, 41, 40, 39]);
        assert_eq!(
            read2_supplementary
                .cigar()
                .iter()
                .collect::<Result<Vec<_>, _>>()
                .unwrap(),
            vec![Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 5)]
        );

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn run_restores_hardclipped_records_from_gzipped_fastq() {
        let dir = temp_dir("run_gz");
        let fastq_path = dir.join("reads.fastq.gz");
        let input_bam_path = dir.join("input.bam");
        let output_bam_path = dir.join("output.bam");

        write_test_fastq_gz(
            &fastq_path,
            &[
                ("read1", "ACGTRYN", "ABCDEFG"),
                ("read2", "TTACCGG", "HIJKLMN"),
            ],
        );

        let primary = make_record(
            Some("read1"),
            Flags::empty(),
            cigar(&[(Kind::Match, 7)]),
            b"ACGTRYN",
            &[32, 33, 34, 35, 36, 37, 38],
        );
        let secondary = make_record(
            Some("read1"),
            Flags::SECONDARY,
            cigar(&[(Kind::HardClip, 1), (Kind::Match, 6)]),
            b"CGTRYN",
            &[33, 34, 35, 36, 37, 38],
        );
        let supplementary = make_record(
            Some("read2"),
            Flags::SUPPLEMENTARY | Flags::REVERSE_COMPLEMENTED,
            cigar(&[(Kind::HardClip, 2), (Kind::Match, 5)]),
            b"CCGGT",
            &[39, 38, 37, 36, 35],
        );

        write_test_bam(&input_bam_path, &[primary, secondary, supplementary]);

        let patched_count = run(
            input_bam_path.to_str().unwrap(),
            fastq_path.to_str().unwrap(),
            output_bam_path.to_str().unwrap(),
        )
        .unwrap();

        assert_eq!(patched_count, 2);

        let (_header, records) = read_bam_records(&output_bam_path);
        assert_eq!(records.len(), 3);

        let read1_secondary = records
            .iter()
            .find(|r| {
                let name: &[u8] = r.name().unwrap().as_ref();
                name == b"read1" && r.flags().is_secondary()
            })
            .unwrap();
        assert_eq!(read1_secondary.sequence().as_ref(), b"ACGTRYN");
        assert_eq!(read1_secondary.quality_scores().as_ref(), &[32, 33, 34, 35, 36, 37, 38]);

        let read2_supplementary = records
            .iter()
            .find(|r| {
                let name: &[u8] = r.name().unwrap().as_ref();
                name == b"read2"
            })
            .unwrap();
        assert_eq!(read2_supplementary.sequence().as_ref(), dna::revcomp(b"TTACCGG"));
        assert_eq!(read2_supplementary.quality_scores().as_ref(), &[45, 44, 43, 42, 41, 40, 39]);

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn run_restores_grch38_backed_fixture_end_to_end() {
        let dir = temp_dir("grch38_fixture");
        let input_bam_path = Path::new("tests/fixtures/test.dev.ONT.withHardClips.usingGRCh38seq.bam");
        let fastq_path = dir.join("test.dev.ONT.withHardClips.usingGRCh38seq.fastq");
        let output_bam_path = dir.join("test.dev.ONT.withHardClips.usingGRCh38seq.patched.bam");

        write_fastq_from_primary_records(input_bam_path, &fastq_path);

        let patched_count = run(
            input_bam_path.to_str().unwrap(),
            fastq_path.to_str().unwrap(),
            output_bam_path.to_str().unwrap(),
        )
        .unwrap();

        assert_eq!(patched_count, 20);

        let (_input_header, input_records) = read_bam_records(input_bam_path);
        let (_output_header, output_records) = read_bam_records(&output_bam_path);
        let (fastq_cache, fastq_count) =
            load_fastq_cache(BufReader::new(fs::File::open(&fastq_path).unwrap())).unwrap();

        assert_eq!(fastq_count, 10);
        assert_eq!(input_records.len(), 30);
        assert_eq!(output_records.len(), 30);

        for (input_record, output_record) in input_records.iter().zip(output_records.iter()) {
            assert_eq!(input_record.name(), output_record.name());
            assert_eq!(input_record.flags(), output_record.flags());

            let input_has_hard_clip = has_hard_clip(input_record.cigar()).unwrap();
            let output_has_hard_clip = has_hard_clip(output_record.cigar()).unwrap();
            assert!(!output_has_hard_clip);

            if input_has_hard_clip {
                let read_name: &[u8] = input_record.name().unwrap().as_ref();
                let (fastq_seq, fastq_qual) = fastq_cache.get(read_name).unwrap();
                let expected_seq = if output_record.flags().is_reverse_complemented() {
                    dna::revcomp(fastq_seq)
                } else {
                    fastq_seq.clone()
                };
                let expected_qual = if output_record.flags().is_reverse_complemented() {
                    let mut qual = fastq_qual.clone();
                    qual.reverse();
                    qual
                } else {
                    fastq_qual.clone()
                };

                assert_eq!(output_record.sequence().as_ref(), expected_seq.as_slice());
                assert_eq!(output_record.quality_scores().as_ref(), expected_qual.as_slice());
                assert_eq!(
                    read_consuming_cigar_len(output_record.cigar()).unwrap(),
                    expected_seq.len()
                );
            } else {
                assert_eq!(output_record.sequence().as_ref(), input_record.sequence().as_ref());
                assert_eq!(
                    output_record.quality_scores().as_ref(),
                    input_record.quality_scores().as_ref()
                );
                assert_eq!(
                    output_record
                        .cigar()
                        .iter()
                        .collect::<Result<Vec<_>, _>>()
                        .unwrap(),
                    input_record
                        .cigar()
                        .iter()
                        .collect::<Result<Vec<_>, _>>()
                        .unwrap()
                );
            }

            assert_eq!(
                output_record.sequence().len(),
                output_record.quality_scores().len()
            );
            assert_eq!(
                output_record.sequence().len(),
                read_consuming_cigar_len(output_record.cigar()).unwrap()
            );
        }

        fs::remove_dir_all(dir).unwrap();
    }
}
