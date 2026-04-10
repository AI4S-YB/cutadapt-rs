/// Routines for printing a report.
///
/// Port of cutadapt/report.py.

use std::collections::HashMap;
use std::fmt;
use std::fmt::Write as FmtWrite;

use serde_json::{json, Value};

use crate::adapters::statistics::{AdapterStatistics, EndStatistics};
use crate::statistics::ReadLengthStatistics;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Map of filter names to human-readable descriptions.
pub const FILTERS: &[(&str, &str)] = &[
    ("too_short", "that were too short"),
    ("too_long", "that were too long"),
    ("too_many_n", "with too many N"),
    ("too_many_expected_errors", "with too many exp. errors"),
    ("casava_filtered", "failed CASAVA filter"),
    ("discard_trimmed", "discarded as trimmed"),
    ("discard_untrimmed", "discarded as untrimmed"),
];

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Safe division: returns 0.0 when the denominator is zero or the numerator is
/// `None`.
pub fn safe_divide(numerator: Option<u64>, denominator: u64) -> f64 {
    match numerator {
        Some(n) if denominator > 0 => n as f64 / denominator as f64,
        _ => 0.0,
    }
}

/// Add two `Option<u64>` values: if either is `None`, return the other.
pub fn add_if_not_none(a: Option<u64>, b: Option<u64>) -> Option<u64> {
    match (a, b) {
        (None, x) | (x, None) => x,
        (Some(x), Some(y)) => Some(x + y),
    }
}

// ---------------------------------------------------------------------------
// Statistics
// ---------------------------------------------------------------------------

/// Main statistics aggregator for a cutadapt run.
///
/// Mirrors `cutadapt.report.Statistics`.
#[derive(Debug)]
pub struct Statistics {
    pub paired: Option<bool>,
    pub n: u64,
    pub total_bp: [u64; 2],
    pub with_adapters: [Option<u64>; 2],
    pub quality_trimmed_bp: [Option<u64>; 2],
    pub poly_a_trimmed_lengths: [Option<HashMap<usize, u64>>; 2],
    pub filtered: HashMap<String, u64>,
    pub reverse_complemented: Option<u64>,
    pub read_length_statistics: ReadLengthStatistics,
    pub adapter_stats: [Vec<Box<dyn AdapterStatistics>>; 2],
    collected: bool,
}

impl Default for Statistics {
    fn default() -> Self {
        Self::new()
    }
}

impl Statistics {
    pub fn new() -> Self {
        Self {
            paired: None,
            n: 0,
            total_bp: [0, 0],
            with_adapters: [None, None],
            quality_trimmed_bp: [None, None],
            poly_a_trimmed_lengths: [None, None],
            filtered: HashMap::new(),
            reverse_complemented: None,
            read_length_statistics: ReadLengthStatistics::new(),
            adapter_stats: [Vec::new(), Vec::new()],
            collected: false,
        }
    }

    // -- Derived properties --------------------------------------------------

    /// Total base pairs across R1 and R2.
    pub fn total(&self) -> u64 {
        self.total_bp[0] + self.total_bp[1]
    }

    /// Total quality-trimmed base pairs across R1 and R2.
    pub fn quality_trimmed(&self) -> Option<u64> {
        add_if_not_none(self.quality_trimmed_bp[0], self.quality_trimmed_bp[1])
    }

    /// Poly-A trimmed base pairs per read end.
    pub fn poly_a_trimmed_bp(&self) -> (Option<u64>, Option<u64>) {
        let trimmed = |i: usize| -> Option<u64> {
            self.poly_a_trimmed_lengths[i].as_ref().map(|lengths| {
                lengths
                    .iter()
                    .map(|(&length, &count)| length as u64 * count)
                    .sum()
            })
        };
        (trimmed(0), trimmed(1))
    }

    /// Total poly-A trimmed base pairs.
    pub fn poly_a_trimmed(&self) -> Option<u64> {
        let (a, b) = self.poly_a_trimmed_bp();
        add_if_not_none(a, b)
    }

    /// Total written base pairs (sum of R1 and R2).
    pub fn total_written_bp(&self) -> u64 {
        let (bp1, bp2) = self.read_length_statistics.written_bp();
        bp1 + bp2
    }

    /// Number of written reads or read pairs.
    pub fn written(&self) -> u64 {
        self.read_length_statistics.written_reads()
    }

    /// Fraction of reads written.
    pub fn written_fraction(&self) -> f64 {
        safe_divide(Some(self.read_length_statistics.written_reads()), self.n)
    }

    /// Per-end fraction of reads with adapters.
    pub fn with_adapters_fraction(&self) -> [f64; 2] {
        [
            safe_divide(self.with_adapters[0], self.n),
            safe_divide(self.with_adapters[1], self.n),
        ]
    }

    /// Fraction of total bp that were quality-trimmed.
    pub fn quality_trimmed_fraction(&self) -> f64 {
        safe_divide(self.quality_trimmed(), self.total())
    }

    /// Written base pairs per read end.
    pub fn written_bp(&self) -> (u64, u64) {
        self.read_length_statistics.written_bp()
    }

    /// Fraction of total bp that were written.
    pub fn total_written_bp_fraction(&self) -> f64 {
        safe_divide(Some(self.total_written_bp()), self.total())
    }

    /// Fraction of reads that were reverse-complemented.
    pub fn reverse_complemented_fraction(&self) -> f64 {
        safe_divide(self.reverse_complemented, self.n)
    }

    /// Fraction of reads filtered by a specific filter name.
    pub fn filtered_fraction(&self, filter_name: &str) -> f64 {
        safe_divide(self.filtered.get(filter_name).copied(), self.n)
    }

    /// Fraction of total bp that were poly-A trimmed.
    pub fn poly_a_trimmed_fraction(&self) -> f64 {
        safe_divide(self.poly_a_trimmed(), self.total())
    }

    // -- Collection ----------------------------------------------------------

    /// Collect statistics from processing results.
    ///
    /// * `n` -- total number of reads / read pairs processed.
    /// * `total_bp1` -- total base pairs in R1.
    /// * `total_bp2` -- total base pairs in R2, or `None` for single-end.
    /// * `filter_counts` -- `(name, count)` pairs from filter steps.
    /// * `read_length_stats` -- aggregated read-length statistics from sinks.
    pub fn collect(
        &mut self,
        n: u64,
        total_bp1: u64,
        total_bp2: Option<u64>,
        filter_counts: &[(&str, u64)],
        read_length_stats: ReadLengthStatistics,
    ) {
        assert!(!self.collected, "Cannot call Statistics.collect more than once");
        self.n = n;
        self.total_bp[0] = total_bp1;
        match total_bp2 {
            Some(bp2) => {
                self.paired = Some(true);
                self.total_bp[1] = bp2;
            }
            None => {
                self.paired = Some(false);
            }
        }
        for &(name, count) in filter_counts {
            *self.filtered.entry(name.to_string()).or_insert(0) += count;
        }
        self.read_length_statistics += read_length_stats;
        self.collected = true;
    }

    // -- Merging (iadd equivalent) -------------------------------------------

    /// Merge another `Statistics` into this one (for combining parallel results).
    pub fn merge(&mut self, other: &Statistics) {
        self.n += other.n;
        self.read_length_statistics += other.read_length_statistics.clone();

        if self.paired.is_none() {
            self.paired = other.paired;
        } else if self.paired != other.paired {
            panic!("Incompatible Statistics: paired is not equal");
        }

        self.reverse_complemented =
            add_if_not_none(self.reverse_complemented, other.reverse_complemented);

        for (filter_name, &count) in &other.filtered {
            *self.filtered.entry(filter_name.clone()).or_insert(0) += count;
        }

        for i in 0..2 {
            self.total_bp[i] += other.total_bp[i];
            self.with_adapters[i] =
                add_if_not_none(self.with_adapters[i], other.with_adapters[i]);
            self.quality_trimmed_bp[i] =
                add_if_not_none(self.quality_trimmed_bp[i], other.quality_trimmed_bp[i]);

            // Merge poly-A trimmed lengths
            match (&mut self.poly_a_trimmed_lengths[i], &other.poly_a_trimmed_lengths[i]) {
                (None, Some(other_map)) => {
                    self.poly_a_trimmed_lengths[i] = Some(other_map.clone());
                }
                (Some(self_map), Some(other_map)) => {
                    for (&length, &count) in other_map {
                        *self_map.entry(length).or_insert(0) += count;
                    }
                }
                _ => {}
            }

            // Merge adapter stats
            if !self.adapter_stats[i].is_empty() && !other.adapter_stats[i].is_empty() {
                assert_eq!(
                    self.adapter_stats[i].len(),
                    other.adapter_stats[i].len(),
                    "Incompatible Statistics objects (adapter_stats length)"
                );
                for j in 0..self.adapter_stats[i].len() {
                    self.adapter_stats[i][j].merge(other.adapter_stats[i][j].as_ref());
                }
            } else if !other.adapter_stats[i].is_empty() && self.adapter_stats[i].is_empty() {
                // Cannot clone trait objects directly, so we note this limitation.
                // In practice the first Statistics in a merge chain should be the
                // one that already has adapter_stats populated.
            }
        }
    }

    // -- JSON output ---------------------------------------------------------

    /// Return a JSON-compatible `serde_json::Value` representation.
    pub fn as_json(&self, gc_content: f64) -> Value {
        let filtered_json: HashMap<&str, Option<u64>> = FILTERS
            .iter()
            .map(|&(name, _)| (name, self.filtered.get(name).copied()))
            .collect();
        let written_bp = self.read_length_statistics.written_bp();
        let is_paired = self.paired == Some(true);
        let (poly_a_bp0, poly_a_bp1) = self.poly_a_trimmed_bp();

        let adapters_r1: Vec<Value> = self.adapter_stats[0]
            .iter()
            .map(|astats| {
                adapter_statistics_as_json(
                    astats.as_ref(),
                    self.n,
                    gc_content,
                    self.reverse_complemented.is_some(),
                )
            })
            .collect();

        let adapters_r2: Option<Vec<Value>> = if is_paired {
            Some(
                self.adapter_stats[1]
                    .iter()
                    .map(|astats| {
                        adapter_statistics_as_json(
                            astats.as_ref(),
                            self.n,
                            gc_content,
                            self.reverse_complemented.is_some(),
                        )
                    })
                    .collect(),
            )
        } else {
            None
        };

        json!({
            "read_counts": {
                "input": self.n,
                "filtered": filtered_json,
                "output": self.read_length_statistics.written_reads(),
                "reverse_complemented": self.reverse_complemented,
                "read1_with_adapter": self.with_adapters[0],
                "read2_with_adapter": if is_paired { self.with_adapters[1] } else { None },
            },
            "basepair_counts": {
                "input": self.total(),
                "input_read1": self.total_bp[0],
                "input_read2": if is_paired { Some(self.total_bp[1]) } else { None },
                "quality_trimmed": self.quality_trimmed(),
                "quality_trimmed_read1": self.quality_trimmed_bp[0],
                "quality_trimmed_read2": self.quality_trimmed_bp[1],
                "poly_a_trimmed": self.poly_a_trimmed(),
                "poly_a_trimmed_read1": poly_a_bp0,
                "poly_a_trimmed_read2": poly_a_bp1,
                "output": self.total_written_bp(),
                "output_read1": written_bp.0,
                "output_read2": if is_paired { Some(written_bp.1) } else { None },
            },
            "adapters_read1": adapters_r1,
            "adapters_read2": adapters_r2,
            "poly_a_trimmed_read1": poly_a_trimmed_as_json(&self.poly_a_trimmed_lengths[0]),
            "poly_a_trimmed_read2": poly_a_trimmed_as_json(&self.poly_a_trimmed_lengths[1]),
        })
    }
}

/// Convert adapter statistics to JSON value.
fn adapter_statistics_as_json(
    adapter_statistics: &dyn AdapterStatistics,
    n: u64,
    gc_content: f64,
    has_reverse_complemented: bool,
) -> Value {
    let (front_end, back_end) = adapter_statistics.end_statistics();
    let mut ends = Vec::new();
    let mut total_trimmed_reads: u64 = 0;

    for end_opt in [front_end, back_end] {
        match end_opt {
            None => ends.push(Value::Null),
            Some(end_statistics) => {
                let lengths = end_statistics.lengths();
                let total: u64 = lengths.values().sum();

                let eranges = if end_statistics.allows_partial_matches {
                    let er = ErrorRanges::new(
                        end_statistics.effective_length as usize,
                        end_statistics.max_error_rate,
                    );
                    Some(json!(er.lengths()))
                } else {
                    None
                };

                let base_stats = AdjacentBaseStatistics::new(&end_statistics.adjacent_bases);

                let trimmed_lengths: Vec<Value> = histogram_rows(end_statistics, n, gc_content)
                    .map(|row| {
                        json!({
                            "len": row.length,
                            "expect": (row.expect * 10.0).round() / 10.0,
                            "counts": row.error_counts,
                        })
                    })
                    .collect();

                ends.push(json!({
                    "type": end_statistics.adapter_type,
                    "sequence": end_statistics.sequence,
                    "error_rate": end_statistics.max_error_rate,
                    "indels": end_statistics.indels,
                    "error_lengths": eranges,
                    "matches": total,
                    "adjacent_bases": base_stats.as_json(),
                    "dominant_adjacent_base": base_stats.warnbase(),
                    "trimmed_lengths": trimmed_lengths,
                }));

                total_trimmed_reads += total;
            }
        }
    }

    let on_reverse_complement = if has_reverse_complemented {
        Some(adapter_statistics.reverse_complemented())
    } else {
        None
    };

    // Determine if this is a linked adapter by checking if both ends have statistics.
    let is_linked = front_end.is_some() && back_end.is_some();

    json!({
        "name": adapter_statistics.name(),
        "total_matches": total_trimmed_reads,
        "on_reverse_complement": on_reverse_complement,
        "linked": is_linked,
        "five_prime_end": ends.first().cloned().unwrap_or(Value::Null),
        "three_prime_end": ends.get(1).cloned().unwrap_or(Value::Null),
    })
}

/// Convert poly-A trimmed lengths to JSON value.
fn poly_a_trimmed_as_json(poly_a: &Option<HashMap<usize, u64>>) -> Value {
    match poly_a {
        None => Value::Null,
        Some(map) => {
            let mut entries: Vec<(usize, u64)> = map.iter().map(|(&l, &c)| (l, c)).collect();
            entries.sort_by_key(|&(l, _)| l);
            let items: Vec<Value> = entries
                .iter()
                .map(|&(length, count)| json!({"len": length, "count": count}))
                .collect();
            Value::Array(items)
        }
    }
}

// ---------------------------------------------------------------------------
// ErrorRanges
// ---------------------------------------------------------------------------

/// Representation of the lengths up to which a number of errors is allowed
/// for partial adapter matches.
///
/// The entry at index i in the returned list is the length up to which
/// i errors are allowed. For example, the list `[9, 19, 23]` describes that
/// - 0 errors are allowed up to length 9
/// - 1 error is allowed up to length 19
/// - 2 errors are allowed up to length 23
///
/// The last number in the list is always the length of the adapter sequence.
#[derive(Debug, Clone)]
pub struct ErrorRanges {
    pub length: usize,
    pub error_rate: f64,
    lengths_cache: Vec<usize>,
}

impl ErrorRanges {
    pub fn new(length: usize, error_rate: f64) -> Self {
        let lengths_cache = Self::compute_lengths(length, error_rate);
        Self {
            length,
            error_rate,
            lengths_cache,
        }
    }

    fn compute_lengths(length: usize, error_rate: f64) -> Vec<usize> {
        let max_errors = (error_rate * length as f64) as usize;
        let mut lengths: Vec<usize> = (1..=max_errors)
            .map(|errors| (errors as f64 / error_rate) as usize - 1)
            .collect();
        if lengths.is_empty() || *lengths.last().unwrap() < length {
            lengths.push(length);
        }
        lengths
    }

    /// Return the list of lengths.
    pub fn lengths(&self) -> &[usize] {
        &self.lengths_cache
    }
}

impl fmt::Display for ErrorRanges {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut prev = 1;
        let n = self.lengths_cache.len();
        for (errors, &r) in self.lengths_cache[..n - 1].iter().enumerate() {
            write!(f, "{}-{} bp: {}; ", prev, r, errors)?;
            prev = r + 1;
        }
        let last = self.lengths_cache[n - 1];
        if prev == last {
            write!(f, "{} bp: {}", prev, n - 1)?;
        } else {
            write!(f, "{}-{} bp: {}", prev, last, n - 1)?;
        }
        Ok(())
    }
}

/// Format the error ranges line for a report.
fn error_ranges_str(end_statistics: &EndStatistics) -> String {
    let length = end_statistics.effective_length as usize;
    let error_rate = end_statistics.max_error_rate;
    if end_statistics.allows_partial_matches {
        format!(
            "No. of allowed errors:\n{}\n",
            ErrorRanges::new(length, error_rate)
        )
    } else {
        format!(
            "No. of allowed errors: {}\n",
            (error_rate * length as f64) as i32
        )
    }
}

// ---------------------------------------------------------------------------
// HistogramRow
// ---------------------------------------------------------------------------

/// One row in the "trimmed lengths" histogram.
#[derive(Debug, Clone)]
pub struct HistogramRow {
    pub length: usize,
    pub count: u64,
    pub expect: f64,
    pub max_err: i32,
    pub error_counts: Vec<u64>,
}

/// Yield histogram rows for an adapter end.
pub fn histogram_rows(
    end_statistics: &EndStatistics,
    n: u64,
    gc_content: f64,
) -> impl Iterator<Item = HistogramRow> + '_ {
    let d = end_statistics.lengths();
    let match_probabilities = end_statistics.random_match_probabilities(gc_content);
    let seq_len = end_statistics.sequence.len();
    let errors = &end_statistics.errors;
    let max_error_rate = end_statistics.max_error_rate;
    let effective_length = end_statistics.effective_length;

    let mut sorted_lengths: Vec<usize> = d.keys().copied().collect();
    sorted_lengths.sort();

    sorted_lengths.into_iter().map(move |length| {
        let prob_idx = std::cmp::min(seq_len, length);
        let expect = n as f64
            * if prob_idx < match_probabilities.len() {
                match_probabilities[prob_idx]
            } else {
                *match_probabilities.last().unwrap_or(&0.0)
            };
        let count = d[&length];
        let error_dict = errors.get(&length);
        let max_errors = error_dict
            .map(|ed| ed.keys().copied().max().unwrap_or(0))
            .unwrap_or(0);
        let error_counts: Vec<u64> = (0..=max_errors)
            .map(|e| {
                error_dict
                    .and_then(|ed| ed.get(&e).copied())
                    .unwrap_or(0)
            })
            .collect();
        HistogramRow {
            length,
            count,
            expect,
            max_err: (max_error_rate * std::cmp::min(length, effective_length as usize) as f64)
                as i32,
            error_counts,
        }
    })
}

/// Format a histogram as a string.
fn histogram(end_statistics: &EndStatistics, n: u64, gc_content: f64) -> String {
    let mut s = String::from("length\tcount\texpect\tmax.err\terror counts\n");
    for row in histogram_rows(end_statistics, n, gc_content) {
        let error_str: Vec<String> = row.error_counts.iter().map(|e| e.to_string()).collect();
        writeln!(
            s,
            "{}\t{}\t{:.1}\t{}\t{}",
            row.length,
            row.count,
            row.expect,
            row.max_err,
            error_str.join(" ")
        )
        .unwrap();
    }
    s.push('\n');
    s
}

// ---------------------------------------------------------------------------
// AdjacentBaseStatistics
// ---------------------------------------------------------------------------

/// Analyzes base composition near adapter matches.
#[derive(Debug, Clone)]
pub struct AdjacentBaseStatistics {
    pub bases: HashMap<String, u64>,
    fractions: Option<Vec<(String, f64)>>,
    warnbase: Option<String>,
}

impl AdjacentBaseStatistics {
    pub fn new(bases: &HashMap<String, u64>) -> Self {
        let total: u64 = bases.values().sum();
        if total == 0 {
            return Self {
                bases: bases.clone(),
                fractions: None,
                warnbase: None,
            };
        }
        let mut fractions = Vec::new();
        let mut warnbase = None;
        for base in &["A", "C", "G", "T", ""] {
            let key = base.to_string();
            let text = if base.is_empty() {
                "none/other".to_string()
            } else {
                base.to_string()
            };
            let count = bases.get(&key).copied().unwrap_or(0);
            let fraction = count as f64 / total as f64;
            fractions.push((text.clone(), fraction));
            if fraction > 0.8 && !base.is_empty() {
                warnbase = Some(text);
            }
        }
        if total < 20 {
            warnbase = None;
        }
        Self {
            bases: bases.clone(),
            fractions: Some(fractions),
            warnbase,
        }
    }

    /// Whether a warning should be shown about dominant adjacent base.
    pub fn should_warn(&self) -> bool {
        self.warnbase.is_some()
    }

    /// The dominant adjacent base, if any.
    pub fn warnbase(&self) -> Option<&str> {
        self.warnbase.as_deref()
    }

    /// Return a JSON representation of the base counts.
    pub fn as_json(&self) -> Option<Value> {
        if self.fractions.is_some() {
            let mut map = serde_json::Map::new();
            for base in &["A", "C", "G", "T", ""] {
                let key = base.to_string();
                let count = self.bases.get(&key).copied().unwrap_or(0);
                map.insert(key, json!(count));
            }
            Some(Value::Object(map))
        } else {
            None
        }
    }
}

impl fmt::Display for AdjacentBaseStatistics {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(ref fractions) = self.fractions {
            writeln!(f, "Bases preceding removed adapters:")?;
            for (text, fraction) in fractions {
                writeln!(f, "  {}: {:.1}%", text, fraction * 100.0)?;
            }
            if self.should_warn() {
                writeln!(f, "WARNING:")?;
                writeln!(
                    f,
                    "    The adapter is preceded by '{}' extremely often.",
                    self.warnbase.as_deref().unwrap()
                )?;
                writeln!(
                    f,
                    "    The provided adapter sequence could be incomplete at its 5' end."
                )?;
                writeln!(f, "    Ignore this warning when trimming primers.")?;
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// full_report
// ---------------------------------------------------------------------------

/// Generate the full text report for a cutadapt run.
///
/// * `stats` -- collected statistics.
/// * `time` -- processing time in seconds.
/// * `gc_content` -- assumed GC content for expected random matches.
pub fn full_report(stats: &Statistics, _time: f64, gc_content: f64) -> String {
    if stats.n == 0 {
        return "No reads processed!".to_string();
    }

    let mut report = String::new();
    let is_paired = stats.paired == Some(true);

    report.push_str("\n=== Summary ===\n\n");

    if is_paired {
        writeln!(report, "Total read pairs processed:      {:>13}", format_number(stats.n)).unwrap();
        for i in 0..2 {
            if let Some(wa) = stats.with_adapters[i] {
                writeln!(
                    report,
                    "  Read {} with adapter:           {:>13} ({:.1}%)",
                    i + 1,
                    format_number(wa),
                    stats.with_adapters_fraction()[i] * 100.0
                )
                .unwrap();
            }
        }
    } else {
        writeln!(report, "Total reads processed:           {:>13}", format_number(stats.n)).unwrap();
        if let Some(wa) = stats.with_adapters[0] {
            writeln!(
                report,
                "Reads with adapters:             {:>13} ({:.1}%)",
                format_number(wa),
                stats.with_adapters_fraction()[0] * 100.0
            )
            .unwrap();
        }
    }

    if stats.reverse_complemented.is_some() {
        writeln!(
            report,
            "Reverse-complemented:            {:>13} ({:.1}%)",
            format_number(stats.reverse_complemented.unwrap()),
            stats.reverse_complemented_fraction() * 100.0
        )
        .unwrap();
    }

    let filter_report = format_filter_report(stats);
    if !filter_report.is_empty() {
        report.push_str("\n== Read fate breakdown ==\n");
        report.push_str(&filter_report);
    }

    let pairs_or_reads = if is_paired { "Pairs" } else { "Reads" };
    writeln!(
        report,
        "{} written (passing filters): {:>13} ({:.1}%)",
        pairs_or_reads,
        format_number(stats.written()),
        stats.written_fraction() * 100.0
    )
    .unwrap();
    report.push('\n');

    writeln!(
        report,
        "Total basepairs processed: {:>13} bp",
        format_number(stats.total())
    )
    .unwrap();
    if is_paired {
        writeln!(report, "  Read 1: {:>13} bp", format_number(stats.total_bp[0])).unwrap();
        writeln!(report, "  Read 2: {:>13} bp", format_number(stats.total_bp[1])).unwrap();
    }

    if let Some(qt) = stats.quality_trimmed() {
        writeln!(
            report,
            "Quality-trimmed:           {:>13} bp ({:.1}%)",
            format_number(qt),
            stats.quality_trimmed_fraction() * 100.0
        )
        .unwrap();
        if is_paired {
            for i in 0..2 {
                if let Some(qtbp) = stats.quality_trimmed_bp[i] {
                    writeln!(report, "  Read {}: {:>13} bp", i + 1, format_number(qtbp)).unwrap();
                }
            }
        }
    }

    if let Some(pat) = stats.poly_a_trimmed() {
        writeln!(
            report,
            "Poly-A-trimmed:            {:>13} bp ({:.1}%)",
            format_number(pat),
            stats.poly_a_trimmed_fraction() * 100.0
        )
        .unwrap();
        if is_paired {
            let (pa0, pa1) = stats.poly_a_trimmed_bp();
            for (i, pa) in [pa0, pa1].iter().enumerate() {
                if let Some(bp) = pa {
                    writeln!(report, "  Read {}: {:>13} bp", i + 1, format_number(*bp)).unwrap();
                }
            }
        }
    }

    writeln!(
        report,
        "Total written (filtered):  {:>13} bp ({:.1}%)",
        format_number(stats.total_written_bp()),
        stats.total_written_bp_fraction() * 100.0
    )
    .unwrap();
    if is_paired {
        let (bp1, bp2) = stats.written_bp();
        writeln!(report, "  Read 1: {:>13} bp", format_number(bp1)).unwrap();
        writeln!(report, "  Read 2: {:>13} bp", format_number(bp2)).unwrap();
    }

    report.push('\n');

    // Per-adapter statistics
    let mut warning = false;
    for which_in_pair in 0..2 {
        for adapter_statistics in &stats.adapter_stats[which_in_pair] {
            let (front_end, back_end) = adapter_statistics.end_statistics();
            let total_front: u64 = front_end
                .map(|es| es.lengths().values().sum())
                .unwrap_or(0);
            let total_back: u64 = back_end
                .map(|es| es.lengths().values().sum())
                .unwrap_or(0);
            let total = total_front + total_back;
            let reverse_complemented = adapter_statistics.reverse_complemented();

            let extra = if is_paired {
                if which_in_pair == 0 {
                    "First read: "
                } else {
                    "Second read: "
                }
            } else {
                ""
            };

            writeln!(
                report,
                "=== {}Adapter {} ===",
                extra,
                adapter_statistics.name()
            )
            .unwrap();
            report.push('\n');

            // Determine adapter type for display
            let is_linked = front_end.is_some() && back_end.is_some();
            let is_front_only = front_end.is_some() && back_end.is_none();
            let is_back_only = front_end.is_none() && back_end.is_some();

            if is_linked {
                let front_es = front_end.unwrap();
                let back_es = back_end.unwrap();
                write!(
                    report,
                    "Sequence: {}...{}; Type: linked; Length: {}+{}; \
                     5' trimmed: {} times; 3' trimmed: {} times",
                    front_es.sequence,
                    back_es.sequence,
                    front_es.sequence.len(),
                    back_es.sequence.len(),
                    total_front,
                    total_back,
                )
                .unwrap();
            } else {
                let end_stats = front_end.or(back_end).unwrap();
                write!(
                    report,
                    "Sequence: {}; Type: {}; Length: {}; Trimmed: {} times",
                    end_stats.sequence,
                    end_stats.adapter_type,
                    end_stats.sequence.len(),
                    total,
                )
                .unwrap();
            }

            if stats.reverse_complemented.is_some() {
                writeln!(
                    report,
                    "; Reverse-complemented: {} times",
                    reverse_complemented
                )
                .unwrap();
            } else {
                report.push('\n');
            }

            if total == 0 {
                report.push('\n');
                continue;
            }

            if is_linked {
                let front_es = front_end.unwrap();
                let back_es = back_end.unwrap();
                report.push('\n');
                // TODO: min_overlap display requires adapter reference
                write!(report, "{}", error_ranges_str(front_es)).unwrap();
                write!(report, "{}", error_ranges_str(back_es)).unwrap();
                writeln!(report, "Overview of removed sequences at 5' end").unwrap();
                write!(report, "{}", histogram(front_es, stats.n, gc_content)).unwrap();
                report.push('\n');
                writeln!(report, "Overview of removed sequences at 3' end").unwrap();
                write!(report, "{}", histogram(back_es, stats.n, gc_content)).unwrap();
            } else if is_front_only {
                let es = front_end.unwrap();
                report.push('\n');
                write!(report, "{}", error_ranges_str(es)).unwrap();
                writeln!(report, "Overview of removed sequences").unwrap();
                write!(report, "{}", histogram(es, stats.n, gc_content)).unwrap();
            } else if is_back_only {
                let es = back_end.unwrap();
                report.push('\n');
                write!(report, "{}", error_ranges_str(es)).unwrap();
                let base_stats = AdjacentBaseStatistics::new(&es.adjacent_bases);
                warning = warning || base_stats.should_warn();
                write!(report, "{}", base_stats).unwrap();
                writeln!(report, "Overview of removed sequences").unwrap();
                write!(report, "{}", histogram(es, stats.n, gc_content)).unwrap();
            }
        }

        // Poly-A report
        if let Some(ref poly_a) = stats.poly_a_trimmed_lengths[which_in_pair] {
            let which = if is_paired {
                Some(which_in_pair)
            } else {
                None
            };
            write!(report, "{}", poly_a_report(poly_a, which)).unwrap();
        }
    }

    if warning {
        writeln!(report, "WARNING:").unwrap();
        writeln!(
            report,
            "    One or more of your adapter sequences may be incomplete."
        )
        .unwrap();
        writeln!(report, "    Please see the detailed output above.").unwrap();
    }

    report.trim_end().to_string()
}

// ---------------------------------------------------------------------------
// minimal_report
// ---------------------------------------------------------------------------

/// Create a minimal tabular report suitable for concatenation.
pub fn minimal_report(stats: &Statistics, _time: f64, _gc_content: f64) -> String {
    let is_paired = stats.paired == Some(true);
    let written_bp = stats.read_length_statistics.written_bp();

    let mut fields: Vec<String> = vec![
        "OK".to_string(),
        stats.n.to_string(),
        stats.total().to_string(),
        stats.filtered.get("too_short").copied().unwrap_or(0).to_string(),
        stats.filtered.get("too_long").copied().unwrap_or(0).to_string(),
        stats.filtered.get("too_many_n").copied().unwrap_or(0).to_string(),
        stats.read_length_statistics.written_reads().to_string(),
        stats.with_adapters[0].unwrap_or(0).to_string(),
        stats.quality_trimmed_bp[0].unwrap_or(0).to_string(),
        written_bp.0.to_string(),
    ];

    if is_paired {
        fields.push(stats.with_adapters[1].unwrap_or(0).to_string());
        fields.push(stats.quality_trimmed_bp[1].unwrap_or(0).to_string());
        fields.push(written_bp.1.to_string());
    }

    // Check for warnings
    let mut has_warning = false;
    for which_in_pair in 0..2 {
        for adapter_statistics in &stats.adapter_stats[which_in_pair] {
            let (_, back_end) = adapter_statistics.end_statistics();
            if let Some(es) = back_end {
                if AdjacentBaseStatistics::new(&es.adjacent_bases).should_warn() {
                    has_warning = true;
                    break;
                }
            }
        }
        if has_warning {
            break;
        }
    }
    if has_warning {
        fields[0] = "WARN".to_string();
    }

    let mut header = vec![
        "status",
        "in_reads",
        "in_bp",
        "too_short",
        "too_long",
        "too_many_n",
        "out_reads",
        "w/adapters",
        "qualtrim_bp",
        "out_bp",
    ];
    if is_paired {
        header.push("w/adapters2");
        header.push("qualtrim2_bp");
        header.push("out2_bp");
    }

    format!("{}\n{}", header.join("\t"), fields.join("\t"))
}

// ---------------------------------------------------------------------------
// Formatting helpers
// ---------------------------------------------------------------------------

/// Format the filter report section.
fn format_filter_report(stats: &Statistics) -> String {
    let is_paired = stats.paired == Some(true);
    let pairs_or_reads = if is_paired { "Pairs" } else { "Reads" };
    let mut report = String::new();
    for &(name, description) in FILTERS {
        if !stats.filtered.contains_key(name) {
            continue;
        }
        let value = stats.filtered[name];
        let fraction = stats.filtered_fraction(name);
        writeln!(
            report,
            "{} {:27} {:>13} ({:.1}%)",
            pairs_or_reads,
            format!("{}:", description),
            format_number(value),
            fraction * 100.0
        )
        .unwrap();
    }
    report
}

/// Format a poly-A report section.
fn poly_a_report(poly_a: &HashMap<usize, u64>, which_in_pair: Option<usize>) -> String {
    let title = match which_in_pair {
        None => "Poly-A",
        Some(0) => "R1 poly-A",
        Some(_) => "R2 poly-A",
    };
    let mut s = format!("=== {} trimmed ===\n\n", title);
    s.push_str("length\tcount\n");
    let mut sorted_lengths: Vec<usize> = poly_a.keys().copied().collect();
    sorted_lengths.sort();
    for length in sorted_lengths {
        writeln!(s, "{}\t{}", length, poly_a[&length]).unwrap();
    }
    s.push('\n');
    s
}

/// Format a u64 with thousands separators (commas).
fn format_number(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::new();
    for (i, ch) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            result.push(',');
        }
        result.push(ch);
    }
    result.chars().rev().collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapters::statistics::{BackAdapterStatistics, EndStatistics, FrontAdapterStatistics};

    #[test]
    fn test_safe_divide() {
        assert_eq!(safe_divide(Some(10), 100), 0.1);
        assert_eq!(safe_divide(None, 100), 0.0);
        assert_eq!(safe_divide(Some(10), 0), 0.0);
    }

    #[test]
    fn test_add_if_not_none() {
        assert_eq!(add_if_not_none(Some(1), Some(2)), Some(3));
        assert_eq!(add_if_not_none(None, Some(2)), Some(2));
        assert_eq!(add_if_not_none(Some(1), None), Some(1));
        assert_eq!(add_if_not_none(None, None), None);
    }

    #[test]
    fn test_error_ranges_short() {
        let er = ErrorRanges::new(8, 0.1);
        assert_eq!(er.lengths(), &[8]);
        assert_eq!(er.to_string(), "1-8 bp: 0");
    }

    #[test]
    fn test_error_ranges_medium() {
        let er = ErrorRanges::new(19, 0.1);
        assert_eq!(er.lengths(), &[9, 19]);
    }

    #[test]
    fn test_error_ranges_20() {
        let er = ErrorRanges::new(20, 0.1);
        assert_eq!(er.lengths(), &[9, 19, 20]);
        assert_eq!(er.to_string(), "1-9 bp: 0; 10-19 bp: 1; 20 bp: 2");
    }

    #[test]
    fn test_error_ranges_21() {
        let er = ErrorRanges::new(21, 0.1);
        assert_eq!(er.lengths(), &[9, 19, 21]);
        assert_eq!(er.to_string(), "1-9 bp: 0; 10-19 bp: 1; 20-21 bp: 2");
    }

    #[test]
    fn test_format_number() {
        assert_eq!(format_number(0), "0");
        assert_eq!(format_number(999), "999");
        assert_eq!(format_number(1000), "1,000");
        assert_eq!(format_number(1234567), "1,234,567");
    }

    #[test]
    fn test_statistics_new() {
        let stats = Statistics::new();
        assert_eq!(stats.n, 0);
        assert_eq!(stats.total(), 0);
        assert!(stats.paired.is_none());
        assert!(stats.quality_trimmed().is_none());
    }

    #[test]
    fn test_statistics_collect_single_end() {
        let mut stats = Statistics::new();
        let mut rls = ReadLengthStatistics::new();
        rls.update(100);
        rls.update(90);
        stats.collect(
            10,
            1000,
            None,
            &[("too_short", 3), ("too_long", 1)],
            rls,
        );
        assert_eq!(stats.n, 10);
        assert_eq!(stats.paired, Some(false));
        assert_eq!(stats.total_bp[0], 1000);
        assert_eq!(stats.total_bp[1], 0);
        assert_eq!(stats.total(), 1000);
        assert_eq!(stats.filtered.get("too_short"), Some(&3));
        assert_eq!(stats.filtered.get("too_long"), Some(&1));
        assert_eq!(stats.written(), 2);
    }

    #[test]
    fn test_statistics_collect_paired_end() {
        let mut stats = Statistics::new();
        let mut rls = ReadLengthStatistics::new();
        rls.update2(100, 80);
        stats.collect(5, 500, Some(400), &[], rls);
        assert_eq!(stats.n, 5);
        assert_eq!(stats.paired, Some(true));
        assert_eq!(stats.total(), 900);
        assert_eq!(stats.written(), 1);
    }

    #[test]
    fn test_statistics_derived_properties() {
        let mut stats = Statistics::new();
        stats.n = 100;
        stats.total_bp = [5000, 4500];
        stats.with_adapters = [Some(30), Some(20)];
        stats.quality_trimmed_bp = [Some(200), Some(150)];
        stats.paired = Some(true);

        assert_eq!(stats.total(), 9500);
        assert_eq!(stats.quality_trimmed(), Some(350));
        let fracs = stats.with_adapters_fraction();
        assert!((fracs[0] - 0.3).abs() < 1e-10);
        assert!((fracs[1] - 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_statistics_merge() {
        let mut s1 = Statistics::new();
        s1.n = 10;
        s1.paired = Some(false);
        s1.total_bp = [500, 0];
        s1.with_adapters = [Some(3), None];
        s1.quality_trimmed_bp = [Some(50), None];
        s1.filtered.insert("too_short".to_string(), 2);

        let mut s2 = Statistics::new();
        s2.n = 20;
        s2.paired = Some(false);
        s2.total_bp = [1000, 0];
        s2.with_adapters = [Some(7), None];
        s2.quality_trimmed_bp = [Some(100), None];
        s2.filtered.insert("too_short".to_string(), 5);
        s2.filtered.insert("too_long".to_string(), 1);

        s1.merge(&s2);
        assert_eq!(s1.n, 30);
        assert_eq!(s1.total_bp[0], 1500);
        assert_eq!(s1.with_adapters[0], Some(10));
        assert_eq!(s1.quality_trimmed_bp[0], Some(150));
        assert_eq!(s1.filtered["too_short"], 7);
        assert_eq!(s1.filtered["too_long"], 1);
    }

    #[test]
    fn test_statistics_merge_poly_a() {
        let mut s1 = Statistics::new();
        s1.paired = Some(false);
        let mut map1 = HashMap::new();
        map1.insert(5, 10);
        map1.insert(10, 3);
        s1.poly_a_trimmed_lengths[0] = Some(map1);

        let mut s2 = Statistics::new();
        s2.paired = Some(false);
        let mut map2 = HashMap::new();
        map2.insert(5, 7);
        map2.insert(15, 2);
        s2.poly_a_trimmed_lengths[0] = Some(map2);

        s1.merge(&s2);
        let merged = s1.poly_a_trimmed_lengths[0].as_ref().unwrap();
        assert_eq!(merged[&5], 17);
        assert_eq!(merged[&10], 3);
        assert_eq!(merged[&15], 2);
    }

    #[test]
    fn test_adjacent_base_statistics_empty() {
        let mut bases = HashMap::new();
        for b in &["A", "C", "G", "T", ""] {
            bases.insert(b.to_string(), 0);
        }
        let abs = AdjacentBaseStatistics::new(&bases);
        assert!(!abs.should_warn());
        assert!(abs.as_json().is_none());
    }

    #[test]
    fn test_adjacent_base_statistics_warn() {
        let mut bases = HashMap::new();
        bases.insert("A".to_string(), 90);
        bases.insert("C".to_string(), 5);
        bases.insert("G".to_string(), 3);
        bases.insert("T".to_string(), 2);
        bases.insert(String::new(), 0);
        let abs = AdjacentBaseStatistics::new(&bases);
        assert!(abs.should_warn());
        assert_eq!(abs.warnbase(), Some("A"));
    }

    #[test]
    fn test_adjacent_base_statistics_no_warn_low_count() {
        let mut bases = HashMap::new();
        bases.insert("A".to_string(), 18);
        bases.insert("C".to_string(), 1);
        bases.insert("G".to_string(), 0);
        bases.insert("T".to_string(), 0);
        bases.insert(String::new(), 0);
        let abs = AdjacentBaseStatistics::new(&bases);
        // total < 20, so no warning
        assert!(!abs.should_warn());
    }

    #[test]
    fn test_full_report_no_reads() {
        let stats = Statistics::new();
        let report = full_report(&stats, 1.0, 0.5);
        assert_eq!(report, "No reads processed!");
    }

    #[test]
    fn test_full_report_single_end() {
        let mut stats = Statistics::new();
        let mut rls = ReadLengthStatistics::new();
        rls.update(100);
        rls.update(95);
        rls.update(80);
        stats.collect(5, 500, None, &[("too_short", 2)], rls);
        let report = full_report(&stats, 1.0, 0.5);
        assert!(report.contains("Total reads processed:"));
        assert!(report.contains("5"));
        assert!(report.contains("too short"));
    }

    #[test]
    fn test_full_report_paired_end() {
        let mut stats = Statistics::new();
        let mut rls = ReadLengthStatistics::new();
        rls.update2(100, 90);
        stats.collect(3, 300, Some(270), &[], rls);
        let report = full_report(&stats, 1.0, 0.5);
        assert!(report.contains("Total read pairs processed:"));
        assert!(report.contains("Read 1:"));
        assert!(report.contains("Read 2:"));
    }

    #[test]
    fn test_minimal_report_single_end() {
        let mut stats = Statistics::new();
        let mut rls = ReadLengthStatistics::new();
        rls.update(100);
        stats.collect(5, 500, None, &[("too_short", 1)], rls);
        let report = minimal_report(&stats, 1.0, 0.5);
        assert!(report.starts_with("status\t"));
        let lines: Vec<&str> = report.lines().collect();
        assert_eq!(lines.len(), 2);
        let fields: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(fields[0], "OK");
        assert_eq!(fields.len(), 10);
    }

    #[test]
    fn test_minimal_report_paired_end() {
        let mut stats = Statistics::new();
        let mut rls = ReadLengthStatistics::new();
        rls.update2(100, 90);
        stats.collect(3, 300, Some(270), &[], rls);
        let report = minimal_report(&stats, 1.0, 0.5);
        let lines: Vec<&str> = report.lines().collect();
        let fields: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(fields.len(), 13); // extra paired fields
    }

    #[test]
    fn test_as_json_basic() {
        let mut stats = Statistics::new();
        let mut rls = ReadLengthStatistics::new();
        rls.update(100);
        stats.collect(5, 500, None, &[("too_short", 2)], rls);
        let json = stats.as_json(0.5);
        assert!(json["read_counts"]["input"].as_u64() == Some(5));
        assert!(json["basepair_counts"]["input"].as_u64() == Some(500));
    }

    #[test]
    fn test_as_json_with_adapters() {
        let mut stats = Statistics::new();
        let rls = ReadLengthStatistics::new();
        stats.collect(100, 10000, None, &[], rls);
        stats.with_adapters[0] = Some(30);

        let es = EndStatistics::new(0.1, "AGATCGGAAG", 10, false, true, "regular_3prime", true, false);
        let adapter_stats = BackAdapterStatistics::new("test_adapter", es);
        stats.adapter_stats[0].push(Box::new(adapter_stats));

        let json = stats.as_json(0.5);
        let adapters = json["adapters_read1"].as_array().unwrap();
        assert_eq!(adapters.len(), 1);
        assert_eq!(adapters[0]["name"], "test_adapter");
    }

    #[test]
    fn test_full_report_with_adapter() {
        let mut stats = Statistics::new();
        let mut rls = ReadLengthStatistics::new();
        rls.update(100);
        rls.update(95);
        stats.collect(10, 1000, None, &[], rls);
        stats.with_adapters[0] = Some(5);

        let mut es = EndStatistics::new(
            0.1,
            "AGATCGGAAG",
            10,
            false,
            true,
            "regular_3prime",
            true,
            false,
        );
        es.add_error(10, 0);
        es.add_error(10, 0);
        es.add_error(8, 0);
        es.add_adjacent_base("A");
        es.add_adjacent_base("A");
        es.add_adjacent_base("C");

        let adapter_stats = BackAdapterStatistics::new("test_adapter", es);
        stats.adapter_stats[0].push(Box::new(adapter_stats));

        let report = full_report(&stats, 1.0, 0.5);
        assert!(report.contains("Adapter test_adapter"));
        assert!(report.contains("AGATCGGAAG"));
        assert!(report.contains("Trimmed: 3 times"));
        assert!(report.contains("Overview of removed sequences"));
    }

    #[test]
    fn test_full_report_with_front_adapter() {
        let mut stats = Statistics::new();
        let mut rls = ReadLengthStatistics::new();
        rls.update(90);
        stats.collect(5, 500, None, &[], rls);
        stats.with_adapters[0] = Some(2);

        let mut es = EndStatistics::new(
            0.1,
            "AGATCGGAAG",
            10,
            false,
            true,
            "regular_5prime",
            true,
            true,
        );
        es.add_error(10, 0);
        es.add_error(8, 1);

        let adapter_stats = FrontAdapterStatistics::new("front_adapter", es);
        stats.adapter_stats[0].push(Box::new(adapter_stats));

        let report = full_report(&stats, 1.0, 0.5);
        assert!(report.contains("Adapter front_adapter"));
        assert!(report.contains("Trimmed: 2 times"));
    }

    #[test]
    fn test_histogram_rows_basic() {
        let mut es = EndStatistics::new(0.1, "AGATCGGAAG", 10, false, true, "test", true, false);
        es.add_error(5, 0);
        es.add_error(5, 0);
        es.add_error(10, 1);

        let rows: Vec<HistogramRow> = histogram_rows(&es, 100, 0.5).collect();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].length, 5);
        assert_eq!(rows[0].count, 2);
        assert_eq!(rows[1].length, 10);
        assert_eq!(rows[1].count, 1);
    }

    #[test]
    fn test_poly_a_trimmed_bp() {
        let mut stats = Statistics::new();
        let mut map = HashMap::new();
        map.insert(5, 10);   // 5 * 10 = 50
        map.insert(10, 3);   // 10 * 3 = 30
        stats.poly_a_trimmed_lengths[0] = Some(map);

        let (bp0, bp1) = stats.poly_a_trimmed_bp();
        assert_eq!(bp0, Some(80));
        assert_eq!(bp1, None);
        assert_eq!(stats.poly_a_trimmed(), Some(80));
    }

    #[test]
    fn test_poly_a_report() {
        let mut map = HashMap::new();
        map.insert(5, 10);
        map.insert(10, 3);
        let report = poly_a_report(&map, None);
        assert!(report.contains("Poly-A trimmed"));
        assert!(report.contains("5\t10"));
        assert!(report.contains("10\t3"));
    }

    #[test]
    fn test_filter_report() {
        let mut stats = Statistics::new();
        stats.n = 100;
        stats.paired = Some(false);
        stats.filtered.insert("too_short".to_string(), 10);
        stats.filtered.insert("too_long".to_string(), 5);
        let report = format_filter_report(&stats);
        assert!(report.contains("too short"));
        assert!(report.contains("too long"));
    }

    #[test]
    fn test_error_ranges_display_single_range() {
        let er = ErrorRanges::new(5, 0.1);
        assert_eq!(er.to_string(), "1-5 bp: 0");
    }
}
