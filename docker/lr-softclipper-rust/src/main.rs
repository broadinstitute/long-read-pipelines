use clap;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

/// Convenience struct for grouping information we extract from a reference
struct RefData {
    // Map of read query names to soft clip lengths in ref records
    soft_clips: HashMap<String, (usize, usize)>,
    // Unmapped reads in ref records
    unmapped_reads: HashSet<String>,
    // Primary alignments from the reference
    alignments: HashMap<String, bam::Record>,
    // Header for the reference bam
    header: bam::Header,
}

/// Takes the `flags` from a bam/sam/cram record and returns true if the read is a primary read, and
/// false otherwise
fn is_primary(flags: u16) -> bool {
    return (flags & (1 << 8)) + (flags & (1 << 11)) == 0;
}

/// Opens the reference bam/sam/cram at `ref_path`, reads through the records in it, and returns a
/// RefData instance containing data from those reads
fn process_reference(ref_path: &Path) -> RefData {
    // Keep track of the soft clips
    let mut soft_clips: HashMap<String, (usize, usize)> = HashMap::new();
    // Keep a list of the unmapped reads in the reference
    let mut unmapped_reads: HashSet<String> = HashSet::new();
    // Create a reader for the ref file
    let mut ref_alignments_file: bam::Reader =
        bam::Reader::from_path(&ref_path).expect("Failed to open reference file for reading.");
    // Keep track of alignments from the reference
    let mut alignments: HashMap<String, bam::Record> = HashMap::new();

    // Loop through records in the reference and build lists of alignments, soft clips, and
    // unmapped reads
    let mut record: bam::Record = bam::Record::new();
    while let Some(read_result) = ref_alignments_file.read(&mut record) {
        // If there was an error reading the record, panic
        if let Err(e) = read_result {
            panic!("Failed to read record from ref, with error {}", e);
        }
        // If this is not a primary read, skip it
        if !is_primary(record.flags()) {
            continue;
        }
        // Parse the query name of the record as a string
        let record_qname: String = match std::str::from_utf8(record.qname()) {
            Ok(qname_str) => String::from(qname_str),
            Err(e) => panic!(
                "Failed to parse record queryname {:?} as string with error {}",
                record.qname(), e
            ),
        };
        // Add the record to our ref alignment list
        alignments.insert(record_qname.clone(), record.clone());
        // If it's unmapped, also add it to our list of unmapped reads
        if record.is_unmapped() {
            unmapped_reads.insert(record_qname.clone());
        }
        // Get the cigar and use it to check for soft clips
        let mut left_clip_len: usize = record.cigar().leading_softclips() as usize;
        let mut right_clip_len: usize = record.cigar().trailing_softclips() as usize;
        // Make sure reads are aligned in the same direction
        if record.is_reverse() {
            let swap = left_clip_len;
            left_clip_len = right_clip_len;
            right_clip_len = swap;
        }
        // Keep track of the soft clips for these reads
        soft_clips.insert(record_qname, (left_clip_len, right_clip_len));
    }

    // Return ref data
    RefData {
        soft_clips,
        unmapped_reads,
        alignments,
        header: bam::Header::from_template(ref_alignments_file.header()),
    }
}

/// Creates 2 bam writers for logging ref conflicts, one for the records from the ref and one for
/// the records from the input file.  Prefixes the file names with `ref_conflicts_prefix` and uses
/// the supplied headers for the files
fn create_writers_for_ref_conflicts(
    ref_conflicts_prefix: &str,
    input_header: &bam::Header,
    ref_header: &bam::Header,
) -> Result<(bam::Writer, bam::Writer), rust_htslib::errors::Error> {
    // Build paths for the two files
    let unassisted_bam_path: PathBuf =
        PathBuf::from(format!("{}_unassisted.bam", ref_conflicts_prefix));
    let ref_bam_path: PathBuf = PathBuf::from(format!("{}_ref.bam", ref_conflicts_prefix));
    // Create and return the writers with the headers from the input and ref files
    let unassisted_bam =
        bam::Writer::from_path(&unassisted_bam_path, input_header, bam::Format::Bam)?;
    let ref_bam = bam::Writer::from_path(&ref_bam_path, ref_header, bam::Format::Bam)?;
    Ok((unassisted_bam, ref_bam))
}

/// Creates and returns a clap App for processing the expected cli arguments
fn configure_cli_args() -> clap::App<'static, 'static> {
    clap::App::new("Soft Clipper")
        .arg(
            clap::Arg::with_name("clipping-threshold")
                .long("clipping-threshold")
                .takes_value(true)
                .value_name("INT")
                .required(true)
                .help("(Required) Minimum soft clip length for read splitting"),
        )
        .arg(
            clap::Arg::with_name("split-read-prefix")
                .long("split-read-prefix")
                .value_name("STRING")
                .takes_value(true)
                .required(true)
                .help(
                    "(Required) Split reads will be prefixed with this arguments. Ex. read_name = read123, \
                    split-read-prefix = x1 would produce read123_x1-l for a split on the left side \
                    and read_123_x1-r for a split on the right side. IMPORTANT - this value should \
                    be different on successive rounds of splitting so that there are no duplicate \
                    naming conflicts.",
                ),
        )
        .arg(
            clap::Arg::with_name("ref")
                .long("ref")
                .value_name("FILE")
                .takes_value(true)
                .help("Assisting reference"),
        )
        .arg(
            clap::Arg::with_name("ref-diff-threshold")
                .long("ref-diff-threshold")
                .value_name("INT")
                .takes_value(true)
                .help("Minimum clipping distance from reference alignment to ignore clipping"),
        )
        .arg(
            clap::Arg::with_name("write-ref-conflicts-prefix")
                .long("write-ref-conflicts-prefix")
                .takes_value(true)
                .value_name("STRING")
                .help(
                    "If set, the tool will only output alignments that are in conflict with the \
                    assisting reference. All lines will be printed in SAM format. The first line \
                    will be the unassisted alignment while the second line is the assisted \
                    alignment from the reference. If you would like to separate the unassisted \
                    from the assisted alignments, save the output to a file and pipe it through \
                    sed to split the file. AN EXAMPLE COMMAND WILL BE ADDED HERE",
                ),
        )
}

fn main() {

    // Set up expected cli arguments
    let app: clap::App = configure_cli_args();
    // Get values of cli arguments, parsing to the expected format and panicking if any required
    // values are missing (which should not actually happen because clap should handle that)
    let cli_args: clap::ArgMatches = app.get_matches();
    let clipping_threshold: usize = cli_args
        .value_of("clipping-threshold")
        .expect("Failed to get value for clipping-threshold from cli.")
        .parse()
        .expect("Failed to parse clipping-threshold as integer.");
    let split_read_prefix: &str = cli_args
        .value_of("split-read-prefix")
        .expect("Failed to get value for split-read-prefix from cli.");
    let reference: Option<PathBuf> = match cli_args.value_of("ref") {
        Some(path_as_str) => Some(PathBuf::from(path_as_str)),
        None => None,
    };
    let ref_diff_threshold: Option<i64> = match cli_args.value_of("ref-diff-threshold") {
        Some(value_as_str) => Some(
            value_as_str
                .parse()
                .expect("Failed to parse ref-diff-threshold as integer."),
        ),
        None => None,
    };
    let write_ref_conflicts_prefix: Option<&str> = cli_args.value_of("write-ref-conflicts-prefix");

    // Create a reader for the input file
    let mut input_reader: bam::Reader =
        //bam::Reader::from_path(&input_file_path).expect("Failed to open input file for reading.");
        bam::Reader::from_stdin().expect("Failed to open input file for reading.");

    // If we have a reference, process the reference
    let (ref_data_maybe, mut unassisted_bam_maybe, mut ref_bam_maybe): (
        Option<RefData>,
        Option<bam::Writer>,
        Option<bam::Writer>,
    ) = match reference {
        Some(ref_file_path) => {
            // Make sure we have a ref diff threshold, since that's necessary for using the
            // reference
            ref_diff_threshold.expect("A value for ref-diff-threshold is required when running with a reference.");
            // Process the records in the reference
            let ref_data = process_reference(&ref_file_path);
            // If the user provided a value for write-ref-conflicts-prefix, open writers for writing
            // conflicts to
            let (unassisted_bam_maybe, ref_bam_maybe): (Option<bam::Writer>, Option<bam::Writer>) =
                match write_ref_conflicts_prefix {
                    Some(prefix) => {
                        match create_writers_for_ref_conflicts(
                            prefix,
                            &bam::Header::from_template(input_reader.header()),
                            &ref_data.header,
                        ) {
                            Ok((unassisted_bam, ref_bam)) => (Some(unassisted_bam), Some(ref_bam)),
                            Err(e) => panic!(
                                "Failed to create writers for ref conflicts with error: {}",
                                e
                            ),
                        }
                    }
                    None => (None, None),
                };
            (Some(ref_data), unassisted_bam_maybe, ref_bam_maybe)
        }
        None => (None, None, None),
    };

    // Loop through the records in our input reader and process them
    let mut record: bam::Record = bam::Record::new();
    while let Some(read_result) = input_reader.read(&mut record) {
        // If there was an error reading the record, panic
        if let Err(e) = read_result {
            panic!("Failed to read record from input file, with error {}", e);
        }
        // If this is not a primary read, skip it
        if !is_primary(record.flags()) {
            continue;
        }
        // Parse the query name of the record as a string
        let record_qname: String = match std::str::from_utf8(record.qname()) {
            Ok(qname_str) => String::from(qname_str),
            Err(e) => panic!(
                "Failed to parse record queryname {:?} as string with error {}",
                record.qname(), e
            ),
        };
        // Get the read sequence
        let record_seq: Vec<u8> = record.seq().as_bytes();
        // Get the cigar and use it to check for soft clips
        let left_clip_len: usize = record.cigar().leading_softclips() as usize;
        let right_clip_len: usize = record.cigar().trailing_softclips() as usize;
        // Initialize the new start and end indices we'll use for the record
        let mut new_primary_start: usize = 0;
        let mut new_primary_end: usize = record.seq_len();
        // Check for soft clips in the reference if one was provided
        let (ref_left_clip_len, ref_right_clip_len): (Option<usize>, Option<usize>) =
            match &ref_data_maybe {
                Some(ref_data) => {
                    // If it's an unmapped read, return 0 for both clip lengths
                    if ref_data.unmapped_reads.contains(&record_qname) {
                        (Some(0), Some(0))
                    } else {
                        // Otherwise, get the left and right soft clip lengths for this read from the
                        // reference
                        match ref_data.soft_clips.get(&record_qname) {
                            Some(ref_clip_lens) => {
                                // If this record is reversed, return the soft clip lengths flipped
                                if record.is_reverse() {
                                    (Some(ref_clip_lens.1), Some(ref_clip_lens.0))
                                } else {
                                    (Some(ref_clip_lens.0), Some(ref_clip_lens.1))
                                }
                            }
                            None => {
                                panic!(
                                    "Reference does not have record with queryname: {}",
                                    record_qname
                                );
                            }
                        }
                    }
                }
                None => (None, None),
            };
        // Check if the length of the left side softclip is more than the clipping threshold
        if left_clip_len >= clipping_threshold {
            // If we have a reference clip and the length difference between it and the left clip
            // for this record exceeds the diff threshold, this is a ref conflict
            if ref_left_clip_len.is_some() && (ref_left_clip_len.unwrap() as i64 - left_clip_len as i64).abs() >= ref_diff_threshold.unwrap() {
                // If we're writing ref conflicts, write this record to the unassisted_bam and the
                // record from the ref to the ref_bam
                if write_ref_conflicts_prefix.is_some() {
                    if let Some(unassisted_bam) = &mut unassisted_bam_maybe{
                        if let Err(e) = unassisted_bam.write(&record){
                            panic!("Encountered the following error ({}) while trying to write record ({:?}) to the ref conflicts unassisted bam", e, record);
                        }
                    }
                    if let Some(ref_bam) = &mut ref_bam_maybe{
                        if let Err(e) = ref_bam.write(ref_data_maybe.as_ref().unwrap().alignments.get(&record_qname).unwrap()) {
                            panic!("Encountered the following error ({}) while trying to write record ({:?}) to the ref conflicts ref bam", e, record);
                        }
                    }
                }
            }
            // If we don't have a reference left clip for this record (or we do and it's not a
            // conflict) print the left clip and set the read to start instead after the clip
            else {
                // Get the subsequence and qual for the left clip
                let sequence: &str = std::str::from_utf8(&record_seq[..left_clip_len]).expect("Failed to parse sequence to string");
                // rust-htslib gives the actual qualities, so we need to add 33 to them to get them
                // in the format they're expected
                let quals: &str = &record.qual()[..left_clip_len].iter().map(|b| (*b + 33) as char).collect::<String>();
                println!("@{}_{}-l\n{}\n+\n{}", record_qname, split_read_prefix, sequence, quals);
                // Set the new start of the read to after the clip
                new_primary_start = left_clip_len as usize;
            }
        }
        // Now we'll do the same thing for right soft clips
        if right_clip_len >= clipping_threshold {
            // If we have a reference clip and the length difference between it and the right clip
            // for this record exceeds the diff threshold, this is a ref conflict
            if ref_right_clip_len.is_some() && (ref_right_clip_len.unwrap() as i64 - right_clip_len as i64).abs() >= ref_diff_threshold.unwrap() {
                // If we're writing ref conflicts, write this record to the unassisted_bam and the
                // record from the ref to the ref_bam
                if write_ref_conflicts_prefix.is_some() {
                    if let Some(unassisted_bam) = &mut unassisted_bam_maybe{
                        if let Err(e) = unassisted_bam.write(&record){
                            panic!("Encountered the following error ({}) while trying to write record ({:?}) to the ref conflicts unassisted bam", e, record);
                        }
                    }
                    if let Some(ref_bam) = &mut ref_bam_maybe{
                        if let Err(e) = ref_bam.write(ref_data_maybe.as_ref().unwrap().alignments.get(&record_qname).unwrap()) {
                            panic!("Encountered the following error ({}) while trying to write record ({:?}) to the ref conflicts ref bam", e, record);
                        }
                    }
                }
            }
            // If we don't have a reference right clip for this record (or we do and it's not a
            // conflict) print the right clip and set the read to start instead after the clip
            else {
                // Get the subsequence and qual for the left clip
                let sequence: &str = std::str::from_utf8(&record_seq[(record.seq_len() - right_clip_len)..]).expect("Failed to parse sequence to string");
                // rust-htslib gives the actual qualities, so we need to add 33 to them to get them
                // in the format they're expected
                let quals: &str = &record.qual()[(record.seq_len() - right_clip_len)..].iter().map(|b| (*b + 33) as char).collect::<String>();
                println!("@{}_{}-r\n{}\n+\n{}", record_qname, split_read_prefix, sequence, quals);
                // Set the new end of the read to before the clip
                new_primary_end = (record.seq_len() - right_clip_len) as usize;
            }
        }
        // Print the read sequence and qual with new start and end
        let sequence: &str = std::str::from_utf8(&record_seq[new_primary_start..new_primary_end]).expect("Failed to parse sequence to string");
        let quals: &str = &record.qual()[new_primary_start..new_primary_end].iter().map(|b| (*b + 33) as char).collect::<String>();
        println!("@{}\n{}\n+\n{}", record_qname, sequence, quals);
    }

}
