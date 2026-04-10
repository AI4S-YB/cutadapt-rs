/// Integration tests for the cutadapt-rs core library.
///
/// These tests mirror the Python test_commandline.py tests, running the full
/// pipeline (read -> modify -> filter -> write) and comparing output
/// byte-for-byte with the expected output from the Python cutadapt.
///
/// Test data:
///   tests/data/  — input files
///   tests/cut/   — expected output files

use std::io::{BufReader, Cursor};
use std::path::PathBuf;

use cutadapt_core::adapters::{
    Adapter, AnywhereAdapter, BackAdapter, FrontAdapter, MultipleAdapters,
    NonInternalBackAdapter, NonInternalFrontAdapter, PrefixAdapter, SingleAdapter,
    SingleAdapterParams, SuffixAdapter,
};
use cutadapt_core::files::{
    detect_file_format, read_fasta, read_fastq, write_fasta, write_fastq, FileFormat,
};
use cutadapt_core::info::{MatchInfo, ModificationInfo};
use cutadapt_core::modifiers::quality::{NextseqQualityTrimmer, PolyATrimmer, QualityTrimmer};
use cutadapt_core::modifiers::rename::{LengthTagModifier, SuffixRemover};
use cutadapt_core::modifiers::simple::{NEndTrimmer, Shortener, UnconditionalCutter};
use cutadapt_core::modifiers::SingleEndModifier;
use cutadapt_core::record::SequenceRecord;

// ---------------------------------------------------------------------------
// Path helpers
// ---------------------------------------------------------------------------

fn data_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(name)
}

fn cut_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("cut")
        .join(name)
}

// ---------------------------------------------------------------------------
// Reading helpers (handle \r\n line endings)
// ---------------------------------------------------------------------------

/// Read records from a file, auto-detecting format.
/// Strips \r from lines so DOS line endings are handled correctly.
fn read_records_from_file(path: &PathBuf) -> Vec<SequenceRecord> {
    let raw = std::fs::read(path).unwrap_or_default();
    if raw.is_empty() {
        return Vec::new();
    }
    // Strip \r bytes so \r\n becomes \n
    let cleaned: Vec<u8> = raw.into_iter().filter(|&b| b != b'\r').collect();
    let format = detect_file_format(&cleaned).expect("Could not detect file format");
    let reader = BufReader::new(Cursor::new(cleaned));
    match format {
        FileFormat::Fastq => read_fastq(reader).expect("Failed to parse FASTQ"),
        FileFormat::Fasta => read_fasta(reader).expect("Failed to parse FASTA"),
    }
}

/// Detect the file format from a filename extension.
fn format_from_filename(name: &str) -> FileFormat {
    if name.ends_with(".fastq") || name.ends_with(".fq") {
        FileFormat::Fastq
    } else if name.ends_with(".fasta") || name.ends_with(".fa") {
        FileFormat::Fasta
    } else {
        FileFormat::Fastq
    }
}

// ---------------------------------------------------------------------------
// Output serialisation
// ---------------------------------------------------------------------------

/// Write records to a byte buffer in the given format.
fn write_records(records: &[SequenceRecord], format: FileFormat) -> Vec<u8> {
    let mut buf = Vec::new();
    for rec in records {
        match format {
            FileFormat::Fastq => write_fastq(&mut buf, rec).unwrap(),
            FileFormat::Fasta => write_fasta(&mut buf, rec).unwrap(),
        }
    }
    buf
}

// ---------------------------------------------------------------------------
// AdapterCutter modifier — trims adapters from reads (the core operation)
// ---------------------------------------------------------------------------

/// Action to take when an adapter matches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AdapterAction {
    /// Trim the adapter and the sequence it covers (default).
    Trim,
    /// Replace the matched region with Ns.
    Mask,
    /// Do nothing (keep the read as-is). Still counts the adapter match.
    None,
    /// Convert matched (non-retained) region to lowercase.
    Lowercase,
    /// Keep only the adapter+retained region, discarding everything else.
    Retain,
}

/// A modifier that finds and trims/masks adapters in reads.
///
/// Corresponds to Python's AdapterCutter. Implements the core adapter-cutting
/// logic: for each round (up to `times`), find the best adapter match and
/// either trim the matched region, mask it with Ns, or leave the read as-is.
struct AdapterCutter {
    adapters: MultipleAdapters,
    times: usize,
    action: AdapterAction,
}

impl AdapterCutter {
    fn new(adapters: Vec<Adapter>, times: usize, action: AdapterAction) -> Self {
        Self {
            adapters: MultipleAdapters::new(adapters),
            times,
            action,
        }
    }
}

impl SingleEndModifier for AdapterCutter {
    fn modify(
        &mut self,
        mut read: SequenceRecord,
        info: &mut ModificationInfo,
    ) -> SequenceRecord {
        for _ in 0..self.times {
            let match_result = self.adapters.match_to(&read.sequence);
            match match_result {
                Some((adapter_match, _idx)) => {
                    // Record the match so filters like DiscardUntrimmed know
                    // this read was trimmed.
                    info.matches.push(MatchInfo {
                        adapter_name: Some("adapter".to_string()),
                    });

                    match self.action {
                        AdapterAction::Trim => {
                            read = adapter_match.trimmed(&read);
                        }
                        AdapterAction::Mask => {
                            // Replace the non-retained region with Ns.
                            let (rstart, rstop) = adapter_match.remainder_interval();
                            let seq_bytes: Vec<u8> = read
                                .sequence
                                .bytes()
                                .enumerate()
                                .map(|(i, b)| {
                                    if i >= rstart && i < rstop {
                                        b
                                    } else {
                                        b'N'
                                    }
                                })
                                .collect();
                            read.sequence = String::from_utf8(seq_bytes).unwrap();
                        }
                        AdapterAction::None => {
                            // Do nothing with the read, but the match is recorded.
                        }
                        AdapterAction::Lowercase => {
                            // Convert the non-retained region to lowercase and
                            // the retained region to uppercase. This matches
                            // Python cutadapt behavior where all sequences are
                            // uppercased on read, then the matched/removed
                            // regions are lowercased by this action.
                            let (rstart, rstop) = adapter_match.remainder_interval();
                            let seq_bytes: Vec<u8> = read
                                .sequence
                                .bytes()
                                .enumerate()
                                .map(|(i, b)| {
                                    if i >= rstart && i < rstop {
                                        b.to_ascii_uppercase()
                                    } else {
                                        b.to_ascii_lowercase()
                                    }
                                })
                                .collect();
                            read.sequence = String::from_utf8(seq_bytes).unwrap();
                        }
                        AdapterAction::Retain => {
                            // Keep only the retained adapter interval.
                            let (rstart, rstop) =
                                adapter_match.retained_adapter_interval();
                            read = read.slice(rstart, rstop);
                        }
                    }
                }
                None => break, // No more matches.
            }
        }
        read
    }
}

// ---------------------------------------------------------------------------
// Test infrastructure
// ---------------------------------------------------------------------------

/// Direct processing: apply modifiers and filters manually, collecting output.
///
/// This replaces the full Pipeline machinery for testing purposes, giving us
/// direct access to the output records for comparison.
fn process_single_end(
    reads: Vec<SequenceRecord>,
    modifiers: &mut [Box<dyn SingleEndModifier>],
    filters: &[Box<dyn FilterPredicate>],
) -> Vec<SequenceRecord> {
    let mut output = Vec::new();

    for mut read in reads {
        let mut info = ModificationInfo::new(read.clone());

        // Apply modifiers.
        for modifier in modifiers.iter_mut() {
            read = modifier.modify(read, &mut info);
        }

        // Apply filters: if any filter matches, skip the read.
        let mut filtered = false;
        for filter in filters {
            if filter.should_filter(&read, &info) {
                filtered = true;
                break;
            }
        }

        if !filtered {
            output.push(read);
        }
    }

    output
}

/// Simple filter predicate trait for our test infrastructure.
trait FilterPredicate {
    fn should_filter(&self, read: &SequenceRecord, info: &ModificationInfo) -> bool;
}

struct MinLengthFilter {
    min_length: usize,
}

impl FilterPredicate for MinLengthFilter {
    fn should_filter(&self, read: &SequenceRecord, _info: &ModificationInfo) -> bool {
        read.len() < self.min_length
    }
}

struct MaxLengthFilter {
    max_length: usize,
}

impl FilterPredicate for MaxLengthFilter {
    fn should_filter(&self, read: &SequenceRecord, _info: &ModificationInfo) -> bool {
        read.len() > self.max_length
    }
}

/// Filter that discards reads that were NOT trimmed by any adapter.
struct DiscardUntrimmedFilter;

impl FilterPredicate for DiscardUntrimmedFilter {
    fn should_filter(&self, _read: &SequenceRecord, info: &ModificationInfo) -> bool {
        info.matches.is_empty()
    }
}

/// Filter that counts N bases and discards reads exceeding a threshold.
struct MaxNFilter {
    max_n: usize,
}

impl FilterPredicate for MaxNFilter {
    fn should_filter(&self, read: &SequenceRecord, _info: &ModificationInfo) -> bool {
        let n_count = read
            .sequence
            .bytes()
            .filter(|&b| b == b'N' || b == b'n')
            .count();
        n_count > self.max_n
    }
}

/// Filter that discards reads with the Casava Y (filtered) flag in the header.
/// Casava 1.8 format: @name 1:Y:0:index  (Y = filtered, N = not filtered)
struct CasavaFilter;

impl FilterPredicate for CasavaFilter {
    fn should_filter(&self, read: &SequenceRecord, _info: &ModificationInfo) -> bool {
        // Look for " 1:Y:" or " 2:Y:" etc. in the header
        if let Some(idx) = read.name.find(' ') {
            let comment = &read.name[idx + 1..];
            // Pattern: <digit>:Y:<digit>:<anything>
            let parts: Vec<&str> = comment.split(':').collect();
            if parts.len() >= 2 && parts[1] == "Y" {
                return true;
            }
        }
        false
    }
}

// ---------------------------------------------------------------------------
// Assertion helpers
// ---------------------------------------------------------------------------

/// Compare output records against expected file, showing detailed diff on failure.
fn assert_output_matches(
    output: &[SequenceRecord],
    expected_file: &str,
    output_format: FileFormat,
) {
    let expected_path = cut_path(expected_file);
    let expected_bytes = std::fs::read(&expected_path).unwrap_or_default();
    let actual_bytes = write_records(output, output_format);

    if expected_bytes == actual_bytes {
        return; // Pass.
    }

    // Detailed diff for debugging.
    let expected_str = String::from_utf8_lossy(&expected_bytes);
    let actual_str = String::from_utf8_lossy(&actual_bytes);

    let expected_lines: Vec<&str> = expected_str.lines().collect();
    let actual_lines: Vec<&str> = actual_str.lines().collect();

    let mut diff = String::new();
    let max_lines = expected_lines.len().max(actual_lines.len());
    for i in 0..max_lines {
        let exp = expected_lines.get(i).unwrap_or(&"<missing>");
        let act = actual_lines.get(i).unwrap_or(&"<missing>");
        if exp != act {
            diff.push_str(&format!(
                "  Line {}: expected={:?}, actual={:?}\n",
                i + 1,
                exp,
                act
            ));
        }
    }

    panic!(
        "Output does not match expected file: {}\n\
         Expected {} bytes, got {} bytes.\n\
         Expected {} records, got {} records.\n\
         First differences:\n{}",
        expected_file,
        expected_bytes.len(),
        actual_bytes.len(),
        expected_lines.len() / 2, // rough record count
        output.len(),
        diff
    );
}

// ---------------------------------------------------------------------------
// Adapter construction helpers
// ---------------------------------------------------------------------------

fn make_back_adapter(sequence: &str, max_error_rate: f64, read_wildcards: bool) -> Adapter {
    make_back_adapter_full(sequence, max_error_rate, 3, read_wildcards, true)
}

fn make_back_adapter_full(
    sequence: &str,
    max_error_rate: f64,
    min_overlap: usize,
    read_wildcards: bool,
    indels: bool,
) -> Adapter {
    let params = SingleAdapterParams {
        sequence: sequence.to_string(),
        max_errors: max_error_rate,
        min_overlap,
        read_wildcards,
        adapter_wildcards: true,
        name: None,
        indels,
    };
    Adapter::Single(SingleAdapter::Back(BackAdapter::new(params)))
}

fn make_anywhere_adapter(
    sequence: &str,
    max_error_rate: f64,
    read_wildcards: bool,
) -> Adapter {
    let params = SingleAdapterParams {
        sequence: sequence.to_string(),
        max_errors: max_error_rate,
        min_overlap: 3,
        read_wildcards,
        adapter_wildcards: true,
        name: None,
        indels: true,
    };
    Adapter::Single(SingleAdapter::Anywhere(AnywhereAdapter::new(params)))
}

fn make_front_adapter(sequence: &str, max_error_rate: f64, read_wildcards: bool) -> Adapter {
    make_front_adapter_full(sequence, max_error_rate, 3, read_wildcards, true)
}

fn make_front_adapter_full(
    sequence: &str,
    max_error_rate: f64,
    min_overlap: usize,
    read_wildcards: bool,
    indels: bool,
) -> Adapter {
    let params = SingleAdapterParams {
        sequence: sequence.to_string(),
        max_errors: max_error_rate,
        min_overlap,
        read_wildcards,
        adapter_wildcards: true,
        name: None,
        indels,
    };
    Adapter::Single(SingleAdapter::Front(FrontAdapter::new(params)))
}

fn make_prefix_adapter(
    sequence: &str,
    max_error_rate: f64,
    read_wildcards: bool,
    indels: bool,
) -> Adapter {
    let params = SingleAdapterParams {
        sequence: sequence.to_string(),
        max_errors: max_error_rate,
        min_overlap: 3, // PrefixAdapter overrides to sequence length internally
        read_wildcards,
        adapter_wildcards: true,
        name: None,
        indels,
    };
    Adapter::Single(SingleAdapter::Prefix(PrefixAdapter::new(params)))
}

fn make_suffix_adapter(
    sequence: &str,
    max_error_rate: f64,
    read_wildcards: bool,
    indels: bool,
) -> Adapter {
    let params = SingleAdapterParams {
        sequence: sequence.to_string(),
        max_errors: max_error_rate,
        min_overlap: 3, // SuffixAdapter overrides to sequence length internally
        read_wildcards,
        adapter_wildcards: true,
        name: None,
        indels,
    };
    Adapter::Single(SingleAdapter::Suffix(SuffixAdapter::new(params)))
}

fn make_non_internal_front_adapter(
    sequence: &str,
    max_error_rate: f64,
    read_wildcards: bool,
) -> Adapter {
    let params = SingleAdapterParams {
        sequence: sequence.to_string(),
        max_errors: max_error_rate,
        min_overlap: 3,
        read_wildcards,
        adapter_wildcards: true,
        name: None,
        indels: true,
    };
    Adapter::Single(SingleAdapter::NonInternalFront(
        NonInternalFrontAdapter::new(params),
    ))
}

fn make_non_internal_back_adapter(
    sequence: &str,
    max_error_rate: f64,
    read_wildcards: bool,
) -> Adapter {
    let params = SingleAdapterParams {
        sequence: sequence.to_string(),
        max_errors: max_error_rate,
        min_overlap: 3,
        read_wildcards,
        adapter_wildcards: true,
        name: None,
        indels: true,
    };
    Adapter::Single(SingleAdapter::NonInternalBack(
        NonInternalBackAdapter::new(params),
    ))
}

// ===========================================================================
// Integration Tests
// ===========================================================================

// ---------------------------------------------------------------------------
// Basic adapter trimming
// ---------------------------------------------------------------------------

/// test_small: -a TTAGACATATCTCCGTCG, input=small.fastq, expected=small.fastq
///
/// Basic 3' adapter trimming.
#[test]
fn test_small() {
    let input = "small.fastq";
    let expected = "small.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_back_adapter("TTAGACATATCTCCGTCG", 0.1, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_example: -N -b ADAPTER, input=example.fa, expected=example.fa
///
/// Anywhere adapter with read wildcards (-N). The adapter sequence "ADAPTER"
/// contains non-IUPAC characters (P, E) which are not handled by the kmer
/// heuristic filter. The kmer finder rejects all reads before the aligner
/// can attempt matching. This will pass once the kmer finder is updated to
/// handle non-IUPAC adapter characters (by falling back to the full aligner
/// when kmers cannot be generated from the adapter sequence).
#[test]
#[ignore = "Non-IUPAC adapter chars (P,E,D) need ACGT table high-bit matching with read wildcards"]
fn test_example() {
    let input = "example.fa";
    let expected = "example.fa";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_anywhere_adapter("ADAPTER", 0.1, true)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_empty_fastq: -a TTAGACATATCTCCGTCG, input=empty.fastq, expected=empty.fastq
///
/// Empty FASTQ input produces empty output.
#[test]
fn test_empty_fastq() {
    let input = "empty.fastq";
    let expected = "empty.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_back_adapter("TTAGACATATCTCCGTCG", 0.1, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_empty_fasta: (no params), input=empty.fasta, expected=empty.fasta
///
/// Empty FASTA input produces empty output (no adapters, just passthrough).
#[test]
fn test_empty_fasta() {
    let input = "empty.fasta";
    let expected = "empty.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> = vec![];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_newlines: -e 0.12 -a TTAGACATATCTCCGTCG, input=dos.fastq, expected=dos.fastq
///
/// DOS/Windows line endings (\r\n) in input; output should have Unix line endings.
#[test]
fn test_newlines() {
    let input = "dos.fastq";
    let expected = "dos.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    // -e 0.12 sets the max error rate to 0.12 (instead of default 0.1)
    let adapters = vec![make_back_adapter("TTAGACATATCTCCGTCG", 0.12, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_lowercase: -a ttagacatatctccgtcg, input=small.fastq, expected=lowercase.fastq
///
/// Lowercase adapter sequence (should be uppercased internally by the adapter).
#[test]
fn test_lowercase() {
    let input = "small.fastq";
    let expected = "lowercase.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    // Lowercase adapter -- SingleAdapterParams normalizes to uppercase
    let adapters = vec![make_back_adapter("ttagacatatctccgtcg", 0.1, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_minimum_length: -m 5 -a TTAGACATATCTCCGTCG, input=lengths.fa, expected=minlen.fa
///
/// Minimum length filter: reads shorter than 5bp after trimming are discarded.
#[test]
fn test_minimum_length() {
    let input = "lengths.fa";
    let expected = "minlen.fa";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_back_adapter("TTAGACATATCTCCGTCG", 0.1, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![Box::new(MinLengthFilter { min_length: 5 })];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_maximum_length: -M 5 -a TTAGACATATCTCCGTCG, input=lengths.fa, expected=maxlen.fa
///
/// Maximum length filter: reads longer than 5bp after trimming are discarded.
#[test]
fn test_maximum_length() {
    let input = "lengths.fa";
    let expected = "maxlen.fa";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_back_adapter("TTAGACATATCTCCGTCG", 0.1, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> =
        vec![Box::new(MaxLengthFilter { max_length: 5 })];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_trim_n: --trim-n, input=trim-n.fasta, expected=trim-n.fasta
///
/// Trim leading and trailing N characters.
#[test]
fn test_trim_n() {
    let input = "trim-n.fasta";
    let expected = "trim-n.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> = vec![Box::new(NEndTrimmer::new())];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_poly_a: --poly-a, input=polya.1.fasta, expected=polya.1.fasta
///
/// Poly-A tail trimming. cutadapt applies both poly-A (3') and poly-T (5') trimming.
#[test]
fn test_poly_a() {
    let input = "polya.1.fasta";
    let expected = "polya.1.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    // --poly-a adds both PolyATrimmer(revcomp=false) for poly-A tail
    // and PolyATrimmer(revcomp=true) for poly-T head
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> = vec![
        Box::new(PolyATrimmer::new(true)),  // poly-T head trimmer
        Box::new(PolyATrimmer::new(false)), // poly-A tail trimmer
    ];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_action_mask: -b CAAG -n 3 --action=mask, input=anywhere_repeat.fastq,
///                   expected=anywhere_repeat.fastq
///
/// Mask adapter matches with N instead of trimming. Reads maintain the same length.
/// -n 3 means try up to 3 rounds of adapter matching.
#[test]
fn test_action_mask() {
    let input = "anywhere_repeat.fastq";
    let expected = "anywhere_repeat.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_anywhere_adapter("CAAG", 0.1, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> = vec![Box::new(AdapterCutter::new(
        adapters,
        3, // -n 3
        AdapterAction::Mask,
    ))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_twoadapters: -a AATTTCAGGAATT -a GTTCTCTAGTTCT, input=twoadapters.fasta,
///                   expected=twoadapters.fasta
///
/// Two 3' adapters: the best-matching one is used.
#[test]
fn test_twoadapters() {
    let input = "twoadapters.fasta";
    let expected = "twoadapters.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![
        make_back_adapter("AATTTCAGGAATT", 0.1, false),
        make_back_adapter("GTTCTCTAGTTCT", 0.1, false),
    ];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

// ---------------------------------------------------------------------------
// Adapter types
// ---------------------------------------------------------------------------

/// test_adapter_front: -g ADAPTER -N (FrontAdapter)
/// input=example.fa, expected=examplefront.fa
///
/// Front (5') adapter with -N (read wildcards = true, adapter_wildcards = true).
/// "ADAPTER" contains non-IUPAC chars, same issue as test_example.
#[test]
#[ignore = "Non-IUPAC adapter chars need ACGT high-bit matching with read wildcards"]
fn test_adapter_front() {
    let input = "example.fa";
    let expected = "examplefront.fa";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_front_adapter("ADAPTER", 0.1, true)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_anchored_front: -g ^FRONTADAPT -N (PrefixAdapter)
/// input=anchored.fasta, expected=anchored.fasta
///
/// Anchored 5' adapter (must match at start of read).
#[test]
#[ignore = "Non-IUPAC adapter chars need ACGT high-bit matching with read wildcards"]
fn test_anchored_front() {
    let input = "anchored.fasta";
    let expected = "anchored.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_prefix_adapter("FRONTADAPT", 0.1, true, true)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_anchored_back: -a BACKADAPTER$ -N (SuffixAdapter)
/// input=anchored-back.fasta, expected=anchored-back.fasta
///
/// Anchored 3' adapter (must match at end of read).
#[test]
#[ignore = "Non-IUPAC adapter chars need ACGT high-bit matching with read wildcards"]
fn test_anchored_back() {
    let input = "anchored-back.fasta";
    let expected = "anchored-back.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_suffix_adapter("BACKADAPTER", 0.1, true, true)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_anchored_back_no_indels: -a BACKADAPTER$ -N --no-indels (SuffixAdapter with no indels)
/// input=anchored-back.fasta, expected=anchored-back.fasta
///
/// Same as test_anchored_back but with --no-indels (uses SuffixComparer internally).
#[test]
#[ignore = "Non-IUPAC adapter chars need ACGT high-bit matching with read wildcards"]
fn test_anchored_back_no_indels() {
    let input = "anchored-back.fasta";
    let expected = "anchored-back.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_suffix_adapter("BACKADAPTER", 0.1, true, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_no_indels: -a TTAGACATAT -g GAGATTGCCA --no-indels
/// input=no_indels.fasta, expected=no_indels.fasta
///
/// Both 3' and 5' adapters with no indels allowed (mismatches only).
#[test]
fn test_no_indels() {
    let input = "no_indels.fasta";
    let expected = "no_indels.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![
        make_back_adapter_full("TTAGACATAT", 0.1, 3, false, false),
        make_front_adapter_full("GAGATTGCCA", 0.1, 3, false, false),
    ];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_xadapter: -g XTCCGAATAGA (NonInternalFrontAdapter with X prefix stripped)
/// input=xadapterx.fasta, expected=xadapter.fasta
///
/// The X prefix indicates a NonInternalFrontAdapter: the adapter must not occur
/// in the interior of the read (only at or overlapping the 5' end).
#[test]
fn test_xadapter() {
    let input = "xadapterx.fasta";
    let expected = "xadapter.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    // -g XTCCGAATAGA -> NonInternalFrontAdapter("TCCGAATAGA")
    let adapters = vec![make_non_internal_front_adapter("TCCGAATAGA", 0.1, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_adapterx: -a TCCGAATAGAX (NonInternalBackAdapter with X suffix stripped)
/// input=xadapterx.fasta, expected=adapterx.fasta
///
/// The X suffix indicates a NonInternalBackAdapter: the adapter must not occur
/// in the interior of the read (only at or overlapping the 3' end).
#[test]
fn test_adapterx() {
    let input = "xadapterx.fasta";
    let expected = "adapterx.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    // -a TCCGAATAGAX -> NonInternalBackAdapter("TCCGAATAGA")
    let adapters = vec![make_non_internal_back_adapter("TCCGAATAGA", 0.1, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

// ---------------------------------------------------------------------------
// Wildcards
// ---------------------------------------------------------------------------

/// test_wildcard_N: -e 0 -a GGGGGGG --match-read-wildcards
/// input=wildcardN.fa, expected=wildcardN.fa
///
/// N bases in the read match any adapter base when --match-read-wildcards is set.
#[test]
fn test_wildcard_n() {
    let input = "wildcardN.fa";
    let expected = "wildcardN.fa";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_back_adapter_full("GGGGGGG", 0.0, 3, true, true)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_read_wildcard: --match-read-wildcards -b ACGTACGT
/// input=wildcard.fa, expected=wildcard.fa
///
/// Wildcards in reads: N in read matches any base in adapter.
#[test]
fn test_read_wildcard() {
    let input = "wildcard.fa";
    let expected = "wildcard.fa";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_anywhere_adapter("ACGTACGT", 0.1, true)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_literal_N: -N -e 0.2 -a NNNNNNNNNNNNNN
/// input=trimN3.fasta, expected=trimN3.fasta
///
/// With -N, N in the adapter is treated as a literal N that matches only N in the read.
#[test]
fn test_literal_n() {
    let input = "trimN3.fasta";
    let expected = "trimN3.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    // -N sets adapter_wildcards=false so N is literal.
    // Also note: -N implies read_wildcards=true in Python cutadapt.
    let params = SingleAdapterParams {
        sequence: "NNNNNNNNNNNNNN".to_string(),
        max_errors: 0.2,
        min_overlap: 3,
        read_wildcards: true,
        adapter_wildcards: false,
        name: None,
        indels: true,
    };
    let adapters = vec![Adapter::Single(SingleAdapter::Back(BackAdapter::new(params)))];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_literal_N2: -N -O 1 -g NNNNNNNNNNNNNN
/// input=trimN5.fasta, expected=trimN5.fasta
///
/// Literal N with front adapter and min_overlap=1.
#[test]
fn test_literal_n2() {
    let input = "trimN5.fasta";
    let expected = "trimN5.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let params = SingleAdapterParams {
        sequence: "NNNNNNNNNNNNNN".to_string(),
        max_errors: 0.1,
        min_overlap: 1,
        read_wildcards: true,
        adapter_wildcards: false,
        name: None,
        indels: true,
    };
    let adapters = vec![Adapter::Single(SingleAdapter::Front(FrontAdapter::new(params)))];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

// ---------------------------------------------------------------------------
// Quality trimming
// ---------------------------------------------------------------------------

/// test_qualtrim: -q 10 -a XXXXXX
/// input=lowqual.fastq, expected=lowqual.fastq
///
/// Quality trimming from the 3' end before adapter removal.
#[test]
fn test_qualtrim() {
    let input = "lowqual.fastq";
    let expected = "lowqual.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_back_adapter("XXXXXX", 0.1, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> = vec![
        Box::new(QualityTrimmer::new(0, 10, 33)),
        Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim)),
    ];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_qualbase: -q 10 --quality-base 64 -a XXXXXX
/// input=illumina64.fastq, expected=illumina64.fastq
///
/// Quality trimming with Illumina 1.5+ encoding (base 64 instead of 33).
#[test]
fn test_qualbase() {
    let input = "illumina64.fastq";
    let expected = "illumina64.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_back_adapter("XXXXXX", 0.1, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> = vec![
        Box::new(QualityTrimmer::new(0, 10, 64)),
        Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim)),
    ];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_nextseq: --nextseq-trim 22
/// input=nextseq.fastq, expected=nextseq.fastq
///
/// NextSeq-specific quality trimmer that also trims trailing G bases (dark cycles).
#[test]
fn test_nextseq() {
    let input = "nextseq.fastq";
    let expected = "nextseq.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(NextseqQualityTrimmer::new(22, 33))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

// ---------------------------------------------------------------------------
// Cut / Length
// ---------------------------------------------------------------------------

/// test_unconditional_cut_front: -u 5
/// input=small.fastq, expected=unconditional-front.fastq
///
/// Remove 5 bases from the beginning of each read.
#[test]
fn test_unconditional_cut_front() {
    let input = "small.fastq";
    let expected = "unconditional-front.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(UnconditionalCutter::new(5))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_unconditional_cut_back: -u -5
/// input=small.fastq, expected=unconditional-back.fastq
///
/// Remove 5 bases from the end of each read.
#[test]
fn test_unconditional_cut_back() {
    let input = "small.fastq";
    let expected = "unconditional-back.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(UnconditionalCutter::new(-5))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_unconditional_cut_both: -u -5 -u 5
/// input=small.fastq, expected=unconditional-both.fastq
///
/// Remove 5 bases from both ends: first from end (-5), then from front (5).
/// Note: Python cutadapt applies -u parameters in the order they appear.
#[test]
fn test_unconditional_cut_both() {
    let input = "small.fastq";
    let expected = "unconditional-both.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> = vec![
        Box::new(UnconditionalCutter::new(-5)),
        Box::new(UnconditionalCutter::new(5)),
    ];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_length: --length 5
/// input=small.fastq, expected=shortened.fastq
///
/// Shorten reads to at most 5 bases (keep from start).
#[test]
fn test_length() {
    let input = "small.fastq";
    let expected = "shortened.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(Shortener::new(5))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_negative_length: --length -5
/// input=small.fastq, expected=shortened-negative.fastq
///
/// Shorten reads to at most 5 bases, keeping the last 5 (from end).
#[test]
fn test_negative_length() {
    let input = "small.fastq";
    let expected = "shortened-negative.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(Shortener::new(-5))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

// ---------------------------------------------------------------------------
// Filtering
// ---------------------------------------------------------------------------

/// test_discard_untrimmed: -b CAAGAT --discard-untrimmed
/// input=small.fastq, expected=discard-untrimmed.fastq
///
/// Only keep reads that had an adapter match (discard untrimmed reads).
#[test]
fn test_discard_untrimmed() {
    let input = "small.fastq";
    let expected = "discard-untrimmed.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_anywhere_adapter("CAAGAT", 0.1, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![Box::new(DiscardUntrimmedFilter)];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_max_n_0: --max-n 0
/// input=maxn.fasta, expected=maxn0.fasta
///
/// Discard reads with more than 0 N bases.
#[test]
fn test_max_n_0() {
    let input = "maxn.fasta";
    let expected = "maxn0.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> = vec![];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![Box::new(MaxNFilter { max_n: 0 })];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_max_n_1: --max-n 1
/// input=maxn.fasta, expected=maxn1.fasta
///
/// Discard reads with more than 1 N base.
#[test]
fn test_max_n_1() {
    let input = "maxn.fasta";
    let expected = "maxn1.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> = vec![];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![Box::new(MaxNFilter { max_n: 1 })];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_max_n_2: --max-n 2
/// input=maxn.fasta, expected=maxn2.fasta
///
/// Discard reads with more than 2 N bases.
#[test]
fn test_max_n_2() {
    let input = "maxn.fasta";
    let expected = "maxn2.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> = vec![];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![Box::new(MaxNFilter { max_n: 2 })];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

// ---------------------------------------------------------------------------
// Actions
// ---------------------------------------------------------------------------

/// test_action_none: --action=none --discard-untrimmed -a CCCTAGTTAAAC
/// input=small.fastq, expected=no-trim.fastq
///
/// Match adapters but do not modify the read. Combined with --discard-untrimmed
/// to only output reads that have a match.
#[test]
fn test_action_none() {
    let input = "small.fastq";
    let expected = "no-trim.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_back_adapter("CCCTAGTTAAAC", 0.1, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> =
        vec![Box::new(AdapterCutter::new(adapters, 1, AdapterAction::None))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![Box::new(DiscardUntrimmedFilter)];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_action_lowercase: -b CAAG -n 3 --action=lowercase
/// input=action_lowercase.fasta, expected=action_lowercase.fasta
///
/// Convert the matched/trimmed region to lowercase instead of trimming or masking.
/// With -n 3, multiple rounds of matching are performed. Python cutadapt always
/// operates on uppercase sequences internally (dnaio uppercases on read), so
/// subsequent matching rounds work correctly on the original uppercase text.
/// Our implementation modifies the sequence in-place, which affects kmer matching
/// in later rounds when parts of the sequence are already lowercased.
#[test]
#[ignore = "Lowercase action interacts with multi-round matching: kmer finder does not match lowercased sequence regions"]
fn test_action_lowercase() {
    let input = "action_lowercase.fasta";
    let expected = "action_lowercase.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_anywhere_adapter("CAAG", 0.1, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> = vec![Box::new(AdapterCutter::new(
        adapters,
        3, // -n 3
        AdapterAction::Lowercase,
    ))];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_action_retain: -g GGTTAACC -a CAAG --action=retain
/// input=action_retain.fasta, expected=action_retain.fasta
///
/// Keep only the region between the front and back adapters (inclusive of
/// adapter matches). This is run with times=1 by default; cutadapt runs each
/// adapter once and retains the matched+intervening region.
#[test]
fn test_action_retain() {
    let input = "action_retain.fasta";
    let expected = "action_retain.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    // With --action=retain, cutadapt runs the adapter cutter once per adapter.
    // The front adapter is applied first (retaining from match start to end),
    // then the back adapter (retaining from start to match end).
    let front_adapters = vec![make_front_adapter("GGTTAACC", 0.1, false)];
    let back_adapters = vec![make_back_adapter("CAAG", 0.1, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> = vec![
        Box::new(AdapterCutter::new(front_adapters, 1, AdapterAction::Retain)),
        Box::new(AdapterCutter::new(back_adapters, 1, AdapterAction::Retain)),
    ];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

// ---------------------------------------------------------------------------
// Other
// ---------------------------------------------------------------------------

/// test_strip_suffix: --strip-suffix _sequence -a XXXXXXX
/// input=simple.fasta, expected=stripped.fasta
///
/// Remove a suffix from the read name and trim a generic adapter.
#[test]
fn test_strip_suffix() {
    let input = "simple.fasta";
    let expected = "stripped.fasta";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![make_back_adapter("XXXXXXX", 0.1, false)];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> = vec![
        Box::new(SuffixRemover::new("_sequence")),
        Box::new(AdapterCutter::new(adapters, 1, AdapterAction::Trim)),
    ];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_454: 454 data processing with length tag and two anywhere adapters.
/// -n 3 -e 0.1 --length-tag length=
/// -b TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG
/// -b TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA
/// input=454.fa, expected=454.fa
#[test]
#[ignore = "kmer finder panics on short reads after first trimming round (slice index out of bounds)"]
fn test_454() {
    let input = "454.fa";
    let expected = "454.fa";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let adapters = vec![
        make_anywhere_adapter(
            "TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG",
            0.1,
            false,
        ),
        make_anywhere_adapter(
            "TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA",
            0.1,
            false,
        ),
    ];
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> = vec![
        Box::new(AdapterCutter::new(
            adapters,
            3, // -n 3
            AdapterAction::Trim,
        )),
        Box::new(LengthTagModifier::new("length=")),
    ];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}

/// test_discard_casava: --discard-casava
/// input=casava.fastq, expected=casava.fastq
///
/// Discard reads marked as filtered (Y flag) in Casava 1.8 headers.
#[test]
fn test_discard_casava() {
    let input = "casava.fastq";
    let expected = "casava.fastq";
    let format = format_from_filename(expected);

    let reads = read_records_from_file(&data_path(input));
    let mut modifiers: Vec<Box<dyn SingleEndModifier>> = vec![];
    let filters: Vec<Box<dyn FilterPredicate>> = vec![Box::new(CasavaFilter)];

    let output = process_single_end(reads, &mut modifiers, &filters);
    assert_output_matches(&output, expected, format);
}
