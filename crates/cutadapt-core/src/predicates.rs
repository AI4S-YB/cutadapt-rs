/// Filtering criteria (predicates).
///
/// Port of cutadapt/predicates.py.
///
/// Each predicate tests whether a read should be filtered out based on some criterion.
/// The `test` method returns `true` if the filtering criterion matches (i.e., the read
/// should be filtered/discarded).

use crate::info::ModificationInfo;
use crate::qualtrim::expected_errors;
use crate::record::SequenceRecord;

/// Default Phred quality encoding base (ASCII 33, i.e. Sanger/Illumina 1.8+).
const DEFAULT_PHRED_BASE: u8 = 33;

/// Trait for all filtering predicates.
pub trait Predicate {
    /// Return `true` if the filtering criterion matches (the read should be filtered).
    fn test(&self, read: &SequenceRecord, info: &ModificationInfo) -> bool;

    /// Return a short snake_case identifier for this predicate, used in reports.
    fn descriptive_identifier(&self) -> &str;
}

/// Select reads that are shorter than the specified minimum length.
#[derive(Debug)]
pub struct TooShort {
    pub minimum_length: usize,
}

impl TooShort {
    pub fn new(minimum_length: usize) -> Self {
        Self { minimum_length }
    }
}

impl Predicate for TooShort {
    fn test(&self, read: &SequenceRecord, _info: &ModificationInfo) -> bool {
        read.len() < self.minimum_length
    }

    fn descriptive_identifier(&self) -> &str {
        "too_short"
    }
}

/// Select reads that are longer than the specified maximum length.
#[derive(Debug)]
pub struct TooLong {
    pub maximum_length: usize,
}

impl TooLong {
    pub fn new(maximum_length: usize) -> Self {
        Self { maximum_length }
    }
}

impl Predicate for TooLong {
    fn test(&self, read: &SequenceRecord, _info: &ModificationInfo) -> bool {
        read.len() > self.maximum_length
    }

    fn descriptive_identifier(&self) -> &str {
        "too_long"
    }
}

/// Select reads that have too many 'N' bases.
///
/// Both a raw count or a proportion (relative to the sequence length) can be used.
/// If `cutoff` is below 1.0, it is treated as a proportion; otherwise it is an absolute count.
#[derive(Debug)]
pub struct TooManyN {
    pub cutoff: f64,
    pub is_proportion: bool,
}

impl TooManyN {
    pub fn new(count: f64) -> Self {
        assert!(count >= 0.0, "count must be non-negative");
        Self {
            is_proportion: count < 1.0,
            cutoff: count,
        }
    }
}

impl Predicate for TooManyN {
    fn test(&self, read: &SequenceRecord, _info: &ModificationInfo) -> bool {
        let n_count = read
            .sequence
            .bytes()
            .filter(|&b| b == b'N' || b == b'n')
            .count();
        if self.is_proportion {
            if read.len() == 0 {
                return false;
            }
            (n_count as f64 / read.len() as f64) > self.cutoff
        } else {
            n_count as f64 > self.cutoff
        }
    }

    fn descriptive_identifier(&self) -> &str {
        "too_many_n"
    }
}

/// Select reads whose expected number of errors, according to quality values,
/// exceeds a threshold.
///
/// The idea comes from usearch's -fastq_maxee parameter
/// (<http://drive5.com/usearch/>).
#[derive(Debug)]
pub struct TooManyExpectedErrors {
    pub max_errors: f64,
}

impl TooManyExpectedErrors {
    pub fn new(max_errors: f64) -> Self {
        Self { max_errors }
    }
}

impl Predicate for TooManyExpectedErrors {
    fn test(&self, read: &SequenceRecord, _info: &ModificationInfo) -> bool {
        let quals = match &read.qualities {
            Some(q) => q.as_bytes(),
            None => return false,
        };
        match expected_errors(quals, DEFAULT_PHRED_BASE) {
            Ok(ee) => ee > self.max_errors,
            Err(_) => false,
        }
    }

    fn descriptive_identifier(&self) -> &str {
        "too_many_expected_errors"
    }
}

/// Select reads that have an average error rate above the threshold.
///
/// This works better than [`TooManyExpectedErrors`] for reads that are expected to
/// have varying lengths, such as for long read sequencing technologies.
#[derive(Debug)]
pub struct TooHighAverageErrorRate {
    pub max_error_rate: f64,
}

impl TooHighAverageErrorRate {
    pub fn new(max_error_rate: f64) -> Self {
        assert!(
            max_error_rate > 0.0 && max_error_rate < 1.0,
            "max_error_rate must be between 0.0 and 1.0 (exclusive), got {}",
            max_error_rate
        );
        Self { max_error_rate }
    }
}

impl Predicate for TooHighAverageErrorRate {
    fn test(&self, read: &SequenceRecord, _info: &ModificationInfo) -> bool {
        let read_length = read.len();
        if read_length == 0 {
            return false;
        }
        let quals = match &read.qualities {
            Some(q) => q.as_bytes(),
            None => return false,
        };
        match expected_errors(quals, DEFAULT_PHRED_BASE) {
            Ok(ee) => (ee / read_length as f64) > self.max_error_rate,
            Err(_) => false,
        }
    }

    fn descriptive_identifier(&self) -> &str {
        "too_high_average_error_rate"
    }
}

/// Select reads that have failed the CASAVA filter according to the read header.
///
/// The headers look like `xxxx x:Y:x:x` (with a `Y`). Reads that pass the filter
/// have an `N` instead of `Y`.
///
/// Reads with unrecognized headers are not selected.
#[derive(Debug)]
pub struct CasavaFiltered;

impl CasavaFiltered {
    pub fn new() -> Self {
        Self
    }
}

impl Default for CasavaFiltered {
    fn default() -> Self {
        Self::new()
    }
}

impl Predicate for CasavaFiltered {
    fn test(&self, read: &SequenceRecord, _info: &ModificationInfo) -> bool {
        // The CASAVA header format: "name index:Y:0:sequence"
        // We look for " " then check characters at positions 1..4 for ":Y:"
        if let Some(space_pos) = read.name.find(' ') {
            let right = &read.name[space_pos + 1..];
            right.len() >= 4 && &right[1..4] == ":Y:"
        } else {
            false
        }
    }

    fn descriptive_identifier(&self) -> &str {
        "casava_filtered"
    }
}

/// Select reads for which no adapter match was found.
#[derive(Debug)]
pub struct IsUntrimmed;

impl IsUntrimmed {
    pub fn new() -> Self {
        Self
    }
}

impl Default for IsUntrimmed {
    fn default() -> Self {
        Self::new()
    }
}

impl Predicate for IsUntrimmed {
    fn test(&self, _read: &SequenceRecord, info: &ModificationInfo) -> bool {
        info.matches.is_empty()
    }

    fn descriptive_identifier(&self) -> &str {
        "discard_untrimmed"
    }
}

/// Select reads for which at least one adapter match was found.
#[derive(Debug)]
pub struct IsTrimmed;

impl IsTrimmed {
    pub fn new() -> Self {
        Self
    }
}

impl Default for IsTrimmed {
    fn default() -> Self {
        Self::new()
    }
}

impl Predicate for IsTrimmed {
    fn test(&self, _read: &SequenceRecord, info: &ModificationInfo) -> bool {
        !info.matches.is_empty()
    }

    fn descriptive_identifier(&self) -> &str {
        "discard_trimmed"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::info::{MatchInfo, ModificationInfo};
    use crate::record::SequenceRecord;

    /// Helper: create a simple read with given sequence and optional qualities.
    fn make_read(name: &str, seq: &str, quals: Option<&str>) -> SequenceRecord {
        SequenceRecord::new(name, seq, quals)
    }

    /// Helper: create a default ModificationInfo for a read.
    fn make_info(read: &SequenceRecord) -> ModificationInfo {
        ModificationInfo::new(read.clone())
    }

    // ── TooShort ──────────────────────────────────────────────────────

    #[test]
    fn test_too_short_matches() {
        let pred = TooShort::new(10);
        let read = make_read("r1", "ACGT", None);
        let info = make_info(&read);
        assert!(pred.test(&read, &info));
    }

    #[test]
    fn test_too_short_no_match() {
        let pred = TooShort::new(4);
        let read = make_read("r1", "ACGT", None);
        let info = make_info(&read);
        assert!(!pred.test(&read, &info));
    }

    #[test]
    fn test_too_short_exact_boundary() {
        let pred = TooShort::new(4);
        let read = make_read("r1", "ACG", None);
        let info = make_info(&read);
        assert!(pred.test(&read, &info)); // 3 < 4
    }

    #[test]
    fn test_too_short_identifier() {
        let pred = TooShort::new(10);
        assert_eq!(pred.descriptive_identifier(), "too_short");
    }

    // ── TooLong ───────────────────────────────────────────────────────

    #[test]
    fn test_too_long_matches() {
        let pred = TooLong::new(3);
        let read = make_read("r1", "ACGT", None);
        let info = make_info(&read);
        assert!(pred.test(&read, &info));
    }

    #[test]
    fn test_too_long_no_match() {
        let pred = TooLong::new(10);
        let read = make_read("r1", "ACGT", None);
        let info = make_info(&read);
        assert!(!pred.test(&read, &info));
    }

    #[test]
    fn test_too_long_exact_boundary() {
        let pred = TooLong::new(4);
        let read = make_read("r1", "ACGT", None);
        let info = make_info(&read);
        assert!(!pred.test(&read, &info)); // 4 is not > 4
    }

    #[test]
    fn test_too_long_identifier() {
        let pred = TooLong::new(10);
        assert_eq!(pred.descriptive_identifier(), "too_long");
    }

    // ── TooManyN ──────────────────────────────────────────────────────

    #[test]
    fn test_too_many_n_absolute_matches() {
        let pred = TooManyN::new(2.0);
        let read = make_read("r1", "ACNNNGT", None);
        let info = make_info(&read);
        assert!(pred.test(&read, &info)); // 3 N's > 2
    }

    #[test]
    fn test_too_many_n_absolute_no_match() {
        let pred = TooManyN::new(5.0);
        let read = make_read("r1", "ACNNNGT", None);
        let info = make_info(&read);
        assert!(!pred.test(&read, &info)); // 3 N's not > 5
    }

    #[test]
    fn test_too_many_n_proportion_matches() {
        let pred = TooManyN::new(0.1);
        let read = make_read("r1", "ANNNNNNNNN", None); // 9 N's out of 10
        let info = make_info(&read);
        assert!(pred.test(&read, &info)); // 0.9 > 0.1
    }

    #[test]
    fn test_too_many_n_proportion_no_match() {
        let pred = TooManyN::new(0.5);
        let read = make_read("r1", "ACGTACGTAN", None); // 1 N out of 10
        let info = make_info(&read);
        assert!(!pred.test(&read, &info)); // 0.1 not > 0.5
    }

    #[test]
    fn test_too_many_n_proportion_empty() {
        let pred = TooManyN::new(0.5);
        let read = make_read("r1", "", None);
        let info = make_info(&read);
        assert!(!pred.test(&read, &info));
    }

    #[test]
    fn test_too_many_n_lowercase() {
        let pred = TooManyN::new(1.0);
        let read = make_read("r1", "ACnnGT", None);
        let info = make_info(&read);
        assert!(pred.test(&read, &info)); // 2 n's > 1
    }

    #[test]
    fn test_too_many_n_identifier() {
        let pred = TooManyN::new(1.0);
        assert_eq!(pred.descriptive_identifier(), "too_many_n");
    }

    // ── TooManyExpectedErrors ─────────────────────────────────────────

    #[test]
    fn test_too_many_expected_errors_matches() {
        let pred = TooManyExpectedErrors::new(0.5);
        // '!' = ASCII 33, phred 0, error rate = 1.0 each
        let read = make_read("r1", "ACGT", Some("!!!!"));
        let info = make_info(&read);
        assert!(pred.test(&read, &info)); // expected errors = 4.0 > 0.5
    }

    #[test]
    fn test_too_many_expected_errors_no_match() {
        let pred = TooManyExpectedErrors::new(10.0);
        // 'I' = ASCII 73, phred 40, error rate ~ 0.0001
        let read = make_read("r1", "ACGT", Some("IIII"));
        let info = make_info(&read);
        assert!(!pred.test(&read, &info));
    }

    #[test]
    fn test_too_many_expected_errors_no_qualities() {
        let pred = TooManyExpectedErrors::new(0.5);
        let read = make_read("r1", "ACGT", None);
        let info = make_info(&read);
        assert!(!pred.test(&read, &info));
    }

    #[test]
    fn test_too_many_expected_errors_identifier() {
        let pred = TooManyExpectedErrors::new(1.0);
        assert_eq!(pred.descriptive_identifier(), "too_many_expected_errors");
    }

    // ── TooHighAverageErrorRate ───────────────────────────────────────

    #[test]
    fn test_too_high_average_error_rate_matches() {
        let pred = TooHighAverageErrorRate::new(0.1);
        // All phred 0 -> error rate 1.0 per base -> average = 1.0 > 0.1
        let read = make_read("r1", "ACGT", Some("!!!!"));
        let info = make_info(&read);
        assert!(pred.test(&read, &info));
    }

    #[test]
    fn test_too_high_average_error_rate_no_match() {
        let pred = TooHighAverageErrorRate::new(0.1);
        // All phred 40 -> average error rate ~ 0.0001 < 0.1
        let read = make_read("r1", "ACGT", Some("IIII"));
        let info = make_info(&read);
        assert!(!pred.test(&read, &info));
    }

    #[test]
    fn test_too_high_average_error_rate_empty() {
        let pred = TooHighAverageErrorRate::new(0.5);
        let read = make_read("r1", "", Some(""));
        let info = make_info(&read);
        assert!(!pred.test(&read, &info));
    }

    #[test]
    fn test_too_high_average_error_rate_no_qualities() {
        let pred = TooHighAverageErrorRate::new(0.1);
        let read = make_read("r1", "ACGT", None);
        let info = make_info(&read);
        assert!(!pred.test(&read, &info));
    }

    #[test]
    #[should_panic(expected = "max_error_rate must be between")]
    fn test_too_high_average_error_rate_invalid_rate() {
        TooHighAverageErrorRate::new(1.5);
    }

    #[test]
    fn test_too_high_average_error_rate_identifier() {
        let pred = TooHighAverageErrorRate::new(0.5);
        assert_eq!(
            pred.descriptive_identifier(),
            "too_high_average_error_rate"
        );
    }

    // ── CasavaFiltered ────────────────────────────────────────────────

    #[test]
    fn test_casava_filtered_matches() {
        let pred = CasavaFiltered::new();
        let read = make_read("read1 1:Y:0:ACGT", "ACGT", None);
        let info = make_info(&read);
        assert!(pred.test(&read, &info));
    }

    #[test]
    fn test_casava_filtered_passes() {
        let pred = CasavaFiltered::new();
        let read = make_read("read1 1:N:0:ACGT", "ACGT", None);
        let info = make_info(&read);
        assert!(!pred.test(&read, &info));
    }

    #[test]
    fn test_casava_filtered_no_space() {
        let pred = CasavaFiltered::new();
        let read = make_read("read1", "ACGT", None);
        let info = make_info(&read);
        assert!(!pred.test(&read, &info));
    }

    #[test]
    fn test_casava_filtered_short_header() {
        let pred = CasavaFiltered::new();
        let read = make_read("read1 X", "ACGT", None);
        let info = make_info(&read);
        assert!(!pred.test(&read, &info)); // right part too short
    }

    #[test]
    fn test_casava_filtered_identifier() {
        let pred = CasavaFiltered::new();
        assert_eq!(pred.descriptive_identifier(), "casava_filtered");
    }

    // ── IsUntrimmed ───────────────────────────────────────────────────

    #[test]
    fn test_is_untrimmed_no_matches() {
        let pred = IsUntrimmed::new();
        let read = make_read("r1", "ACGT", None);
        let info = make_info(&read); // matches is empty
        assert!(pred.test(&read, &info));
    }

    #[test]
    fn test_is_untrimmed_has_matches() {
        let pred = IsUntrimmed::new();
        let read = make_read("r1", "ACGT", None);
        let mut info = make_info(&read);
        info.matches.push(MatchInfo { adapter_name: None });
        assert!(!pred.test(&read, &info));
    }

    #[test]
    fn test_is_untrimmed_identifier() {
        let pred = IsUntrimmed::new();
        assert_eq!(pred.descriptive_identifier(), "discard_untrimmed");
    }

    // ── IsTrimmed ─────────────────────────────────────────────────────

    #[test]
    fn test_is_trimmed_has_matches() {
        let pred = IsTrimmed::new();
        let read = make_read("r1", "ACGT", None);
        let mut info = make_info(&read);
        info.matches.push(MatchInfo { adapter_name: None });
        assert!(pred.test(&read, &info));
    }

    #[test]
    fn test_is_trimmed_no_matches() {
        let pred = IsTrimmed::new();
        let read = make_read("r1", "ACGT", None);
        let info = make_info(&read);
        assert!(!pred.test(&read, &info));
    }

    #[test]
    fn test_is_trimmed_identifier() {
        let pred = IsTrimmed::new();
        assert_eq!(pred.descriptive_identifier(), "discard_trimmed");
    }

    // ── Trait object usage ────────────────────────────────────────────

    #[test]
    fn test_predicates_as_trait_objects() {
        let predicates: Vec<Box<dyn Predicate>> = vec![
            Box::new(TooShort::new(10)),
            Box::new(TooLong::new(100)),
            Box::new(TooManyN::new(5.0)),
            Box::new(TooManyExpectedErrors::new(2.0)),
            Box::new(TooHighAverageErrorRate::new(0.5)),
            Box::new(CasavaFiltered::new()),
            Box::new(IsUntrimmed::new()),
            Box::new(IsTrimmed::new()),
        ];
        let read = make_read("r1", "ACGT", Some("IIII"));
        let info = make_info(&read);
        // Just verify they can all be called through the trait interface
        for pred in &predicates {
            let _ = pred.test(&read, &info);
            let _ = pred.descriptive_identifier();
        }
        assert_eq!(predicates.len(), 8);
    }
}
