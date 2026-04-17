# Testing

Run the Rust test suite with:

```bash
cargo test
```

## Current coverage

The unit and end-to-end tests in `src/main.rs` currently cover:

- `read_consuming_cigar_len` counting only read-consuming CIGAR ops.
- Hard-clip detection on CIGARs with and without `H`.
- Conversion of `H` to `S` while preserving all other CIGAR ops.
- Record classification as `primary`, `secondary`, or `supplementary`.
- FASTQ cache loading, including duplicate read-name overwrite behavior.
- No-op behavior for records without hard clips.
- Forward-strand restoration of sequence, qualities, and CIGAR.
- Reverse-strand restoration of sequence and reversed qualities.
- IUPAC-aware reverse-complement restoration via `bio::alphabets::dna::revcomp`.
- Failure on hard-clipped records missing a read name.
- Failure on hard-clipped records missing a matching FASTQ entry.
- Failure on FASTQ sequence/quality length mismatch.
- Failure when the restored CIGAR read span does not match the restored read length.
- File-based end-to-end restoration through BAM input, FASTQ cache loading, BAM output, and reread verification.
- End-to-end restoration against the committed `tests/fixtures/test.dev.ONT.withHardClips.usingGRCh38seq.bam` fixture, with FASTQ reconstructed from its primary alignments during the test.

## Test intent

The suite is designed to lock down the critical invariants of the program:

- Hard-clipped records are the only records patched.
- Restored `SEQ`, restored `QUAL`, FASTQ sequence length, FASTQ quality length, and read-consuming CIGAR length must all agree.
- Reverse-strand records must use reverse-complemented sequence and reversed qualities.
- Secondary and supplementary records are treated the same as any other hard-clipped record except for classification in diagnostics.

## Layout

Relevant test-related files are organized as follows:

- `src/main.rs`: main binary and unit tests for internal restoration helpers.
- `src/bin/check_fixture_consistency.rs`: manual helper binary for validating fixture BAM consistency.
- `tests/fixtures/`: source BAM fixtures used for manual and automated testing.
- `tests/generated/`: derived artifacts produced from fixtures, such as FASTQ output and patched BAM output.
- `scripts/`: Python generators for regenerating BAM fixtures.
