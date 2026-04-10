/// Steps of the read output pipeline.
///
/// Port of cutadapt/steps.py.
///
/// After all read modifications have been done, a read is written to at
/// most one output file. For this, a pipeline represented as a list of "steps"
/// (SingleEndSteps or PairedEndSteps) is used. Each pipeline step can consume
/// (discard) a read or pass it on to the next step.

use std::collections::HashMap;

use crate::info::ModificationInfo;
use crate::predicates::Predicate;
use crate::record::SequenceRecord;
use crate::statistics::ReadLengthStatistics;

/// A pair of sequence records (R1, R2).
pub type RecordPair = (SequenceRecord, SequenceRecord);

// ---------------------------------------------------------------------------
// Traits
// ---------------------------------------------------------------------------

/// A single-end pipeline step that can consume or pass through a read.
pub trait SingleEndStep {
    /// Process a single read. Return the (possibly modified) read, or `None`
    /// to indicate that the read has been consumed and should not be passed to
    /// subsequent steps.
    fn process(
        &mut self,
        read: SequenceRecord,
        info: &ModificationInfo,
    ) -> Option<SequenceRecord>;
}

/// A paired-end pipeline step that can consume or pass through a read pair.
pub trait PairedEndStep {
    /// Process a read pair. Return the (possibly modified) pair, or `None`
    /// to indicate that the pair has been consumed.
    fn process_pair(
        &mut self,
        read1: SequenceRecord,
        read2: SequenceRecord,
        info1: &ModificationInfo,
        info2: &ModificationInfo,
    ) -> Option<RecordPair>;
}

/// Implemented by final steps (sinks) that track read length statistics.
pub trait HasStatistics {
    fn get_statistics(&self) -> &ReadLengthStatistics;
}

/// Implemented by filtering steps that track how many reads were filtered.
pub trait HasFilterStatistics {
    /// Return the number of filtered reads or read pairs.
    fn filtered(&self) -> u64;

    /// Return a short identifier for this filter, used in reports.
    fn descriptive_identifier(&self) -> &str;
}

// ---------------------------------------------------------------------------
// PairFilterMode
// ---------------------------------------------------------------------------

/// Controls how a paired-end filter decides whether to discard a pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairFilterMode {
    /// Discard the pair if **either** read matches.
    Any,
    /// Discard the pair only if **both** reads match.
    Both,
    /// Discard the pair based on the first read only.
    First,
}

// ---------------------------------------------------------------------------
// SingleEndFilter
// ---------------------------------------------------------------------------

/// A pipeline step that filters reads according to a predicate.
///
/// Optionally collects filtered reads into a separate buffer (standing in for
/// a file writer that will be wired up in Phase 7).
pub struct SingleEndFilter {
    predicate: Box<dyn Predicate>,
    filtered_reads: Option<Vec<SequenceRecord>>,
    filtered_count: u64,
}

impl SingleEndFilter {
    /// Create a new filter. If `collect_filtered` is `true`, filtered reads
    /// are collected into an internal buffer retrievable via
    /// [`Self::take_filtered_reads`].
    pub fn new(predicate: Box<dyn Predicate>, collect_filtered: bool) -> Self {
        Self {
            predicate,
            filtered_reads: if collect_filtered {
                Some(Vec::new())
            } else {
                None
            },
            filtered_count: 0,
        }
    }

    /// Take all collected filtered reads out of the buffer.
    pub fn take_filtered_reads(&mut self) -> Vec<SequenceRecord> {
        self.filtered_reads
            .as_mut()
            .map(std::mem::take)
            .unwrap_or_default()
    }
}

impl SingleEndStep for SingleEndFilter {
    fn process(
        &mut self,
        read: SequenceRecord,
        info: &ModificationInfo,
    ) -> Option<SequenceRecord> {
        if self.predicate.test(&read, info) {
            self.filtered_count += 1;
            if let Some(ref mut buf) = self.filtered_reads {
                buf.push(read);
            }
            None
        } else {
            Some(read)
        }
    }
}

impl HasFilterStatistics for SingleEndFilter {
    fn filtered(&self) -> u64 {
        self.filtered_count
    }

    fn descriptive_identifier(&self) -> &str {
        self.predicate.descriptive_identifier()
    }
}

// ---------------------------------------------------------------------------
// PairedEndFilter
// ---------------------------------------------------------------------------

/// A pipeline step that filters paired-end reads.
///
/// Different filtering styles are supported via [`PairFilterMode`].
pub struct PairedEndFilter {
    predicate1: Option<Box<dyn Predicate>>,
    predicate2: Option<Box<dyn Predicate>>,
    pair_filter_mode: PairFilterMode,
    filtered_pairs: Option<Vec<RecordPair>>,
    filtered_count: u64,
}

impl PairedEndFilter {
    pub fn new(
        predicate1: Option<Box<dyn Predicate>>,
        predicate2: Option<Box<dyn Predicate>>,
        pair_filter_mode: PairFilterMode,
        collect_filtered: bool,
    ) -> Self {
        Self {
            predicate1,
            predicate2,
            pair_filter_mode,
            filtered_pairs: if collect_filtered {
                Some(Vec::new())
            } else {
                None
            },
            filtered_count: 0,
        }
    }

    /// Take all collected filtered pairs out of the buffer.
    pub fn take_filtered_pairs(&mut self) -> Vec<RecordPair> {
        self.filtered_pairs
            .as_mut()
            .map(std::mem::take)
            .unwrap_or_default()
    }

    fn is_filtered(
        &self,
        read1: &SequenceRecord,
        read2: &SequenceRecord,
        info1: &ModificationInfo,
        info2: &ModificationInfo,
    ) -> bool {
        let test1 = self
            .predicate1
            .as_ref()
            .map(|p| p.test(read1, info1))
            .unwrap_or(false);
        let test2 = self
            .predicate2
            .as_ref()
            .map(|p| p.test(read2, info2))
            .unwrap_or(false);

        // If only one predicate is provided, use only that one.
        if self.predicate2.is_none() {
            return test1;
        }
        if self.predicate1.is_none() {
            return test2;
        }

        match self.pair_filter_mode {
            PairFilterMode::Any => test1 || test2,
            PairFilterMode::Both => test1 && test2,
            PairFilterMode::First => test1,
        }
    }
}

impl PairedEndStep for PairedEndFilter {
    fn process_pair(
        &mut self,
        read1: SequenceRecord,
        read2: SequenceRecord,
        info1: &ModificationInfo,
        info2: &ModificationInfo,
    ) -> Option<RecordPair> {
        if self.is_filtered(&read1, &read2, info1, info2) {
            self.filtered_count += 1;
            if let Some(ref mut buf) = self.filtered_pairs {
                buf.push((read1, read2));
            }
            None
        } else {
            Some((read1, read2))
        }
    }
}

impl HasFilterStatistics for PairedEndFilter {
    fn filtered(&self) -> u64 {
        self.filtered_count
    }

    fn descriptive_identifier(&self) -> &str {
        if let Some(ref p) = self.predicate1 {
            p.descriptive_identifier()
        } else if let Some(ref p) = self.predicate2 {
            p.descriptive_identifier()
        } else {
            "unknown"
        }
    }
}

// ---------------------------------------------------------------------------
// SingleEndSink
// ---------------------------------------------------------------------------

/// Final pipeline step for single-end reads.
///
/// Collects reads into an internal buffer and tracks length statistics.
/// (Real file I/O will replace the buffer in Phase 7.)
pub struct SingleEndSink {
    pub reads: Vec<SequenceRecord>,
    statistics: ReadLengthStatistics,
}

impl SingleEndSink {
    pub fn new() -> Self {
        Self {
            reads: Vec::new(),
            statistics: ReadLengthStatistics::new(),
        }
    }
}

impl Default for SingleEndSink {
    fn default() -> Self {
        Self::new()
    }
}

impl SingleEndStep for SingleEndSink {
    fn process(
        &mut self,
        read: SequenceRecord,
        _info: &ModificationInfo,
    ) -> Option<SequenceRecord> {
        self.statistics.update(read.len());
        self.reads.push(read);
        None
    }
}

impl HasStatistics for SingleEndSink {
    fn get_statistics(&self) -> &ReadLengthStatistics {
        &self.statistics
    }
}

// ---------------------------------------------------------------------------
// PairedEndSink
// ---------------------------------------------------------------------------

/// Final pipeline step for paired-end reads.
///
/// Collects read pairs and tracks length statistics.
pub struct PairedEndSink {
    pub pairs: Vec<RecordPair>,
    statistics: ReadLengthStatistics,
}

impl PairedEndSink {
    pub fn new() -> Self {
        Self {
            pairs: Vec::new(),
            statistics: ReadLengthStatistics::new(),
        }
    }
}

impl Default for PairedEndSink {
    fn default() -> Self {
        Self::new()
    }
}

impl PairedEndStep for PairedEndSink {
    fn process_pair(
        &mut self,
        read1: SequenceRecord,
        read2: SequenceRecord,
        _info1: &ModificationInfo,
        _info2: &ModificationInfo,
    ) -> Option<RecordPair> {
        self.statistics.update2(read1.len(), read2.len());
        self.pairs.push((read1, read2));
        None
    }
}

impl HasStatistics for PairedEndSink {
    fn get_statistics(&self) -> &ReadLengthStatistics {
        &self.statistics
    }
}

// ---------------------------------------------------------------------------
// Demultiplexer
// ---------------------------------------------------------------------------

/// Routes single-end reads to different output bins based on the adapter name
/// from the most recent match.
///
/// Reads without a match go to the "untrimmed" bin (if not discarding untrimmed).
pub struct Demultiplexer {
    /// Per-adapter read bins, keyed by adapter name.
    pub bins: HashMap<String, Vec<SequenceRecord>>,
    /// Reads that did not match any adapter.
    pub untrimmed: Option<Vec<SequenceRecord>>,
    statistics: ReadLengthStatistics,
    filtered_count: u64,
}

impl Demultiplexer {
    /// Create a new demultiplexer.
    ///
    /// * `adapter_names` -- the set of expected adapter names.
    /// * `discard_untrimmed` -- if `true`, reads without a match are silently
    ///   discarded; otherwise they are collected in the `untrimmed` bin.
    pub fn new(adapter_names: &[String], discard_untrimmed: bool) -> Self {
        let mut bins = HashMap::new();
        for name in adapter_names {
            bins.insert(name.clone(), Vec::new());
        }
        Self {
            bins,
            untrimmed: if discard_untrimmed {
                None
            } else {
                Some(Vec::new())
            },
            statistics: ReadLengthStatistics::new(),
            filtered_count: 0,
        }
    }
}

impl SingleEndStep for Demultiplexer {
    fn process(
        &mut self,
        read: SequenceRecord,
        info: &ModificationInfo,
    ) -> Option<SequenceRecord> {
        if !info.matches.is_empty() {
            if let Some(ref name) = info.matches.last().and_then(|m| m.adapter_name.as_ref()) {
                self.statistics.update(read.len());
                if let Some(bin) = self.bins.get_mut(name.as_str()) {
                    bin.push(read);
                }
            }
        } else if let Some(ref mut untrimmed) = self.untrimmed {
            self.statistics.update(read.len());
            untrimmed.push(read);
        } else {
            self.filtered_count += 1;
        }
        None
    }
}

impl HasStatistics for Demultiplexer {
    fn get_statistics(&self) -> &ReadLengthStatistics {
        &self.statistics
    }
}

impl HasFilterStatistics for Demultiplexer {
    fn filtered(&self) -> u64 {
        self.filtered_count
    }

    fn descriptive_identifier(&self) -> &str {
        "discard_untrimmed"
    }
}

// ---------------------------------------------------------------------------
// PairedDemultiplexer
// ---------------------------------------------------------------------------

/// Routes paired-end reads based on the adapter name from the R1 match.
pub struct PairedDemultiplexer {
    /// Per-adapter pair bins, keyed by adapter name.
    pub bins: HashMap<String, Vec<RecordPair>>,
    /// Pairs where R1 did not match any adapter.
    pub untrimmed: Option<Vec<RecordPair>>,
    statistics: ReadLengthStatistics,
    filtered_count: u64,
}

impl PairedDemultiplexer {
    pub fn new(adapter_names: &[String], discard_untrimmed: bool) -> Self {
        let mut bins = HashMap::new();
        for name in adapter_names {
            bins.insert(name.clone(), Vec::new());
        }
        Self {
            bins,
            untrimmed: if discard_untrimmed {
                None
            } else {
                Some(Vec::new())
            },
            statistics: ReadLengthStatistics::new(),
            filtered_count: 0,
        }
    }
}

impl PairedEndStep for PairedDemultiplexer {
    fn process_pair(
        &mut self,
        read1: SequenceRecord,
        read2: SequenceRecord,
        info1: &ModificationInfo,
        _info2: &ModificationInfo,
    ) -> Option<RecordPair> {
        if !info1.matches.is_empty() {
            if let Some(ref name) =
                info1.matches.last().and_then(|m| m.adapter_name.as_ref())
            {
                self.statistics.update2(read1.len(), read2.len());
                if let Some(bin) = self.bins.get_mut(name.as_str()) {
                    bin.push((read1, read2));
                }
            }
        } else if let Some(ref mut untrimmed) = self.untrimmed {
            self.statistics.update2(read1.len(), read2.len());
            untrimmed.push((read1, read2));
        } else {
            self.filtered_count += 1;
        }
        None
    }
}

impl HasStatistics for PairedDemultiplexer {
    fn get_statistics(&self) -> &ReadLengthStatistics {
        &self.statistics
    }
}

impl HasFilterStatistics for PairedDemultiplexer {
    fn filtered(&self) -> u64 {
        self.filtered_count
    }

    fn descriptive_identifier(&self) -> &str {
        "discard_untrimmed"
    }
}

// ---------------------------------------------------------------------------
// CombinatorialDemultiplexer
// ---------------------------------------------------------------------------

/// Routes paired-end reads based on the adapter names from both R1 and R2.
///
/// The key is `(Option<adapter_name_r1>, Option<adapter_name_r2>)`.
pub struct CombinatorialDemultiplexer {
    /// Per-(name1, name2) pair bins.
    pub bins: HashMap<(Option<String>, Option<String>), Vec<RecordPair>>,
    statistics: ReadLengthStatistics,
}

impl CombinatorialDemultiplexer {
    /// Create a new combinatorial demultiplexer.
    ///
    /// * `adapter_names1` -- adapter names expected on R1.
    /// * `adapter_names2` -- adapter names expected on R2.
    /// * `discard_untrimmed` -- if `true`, pairs where one or both reads are
    ///   untrimmed are silently discarded; otherwise bins for `None` keys are
    ///   created.
    pub fn new(
        adapter_names1: &[String],
        adapter_names2: &[String],
        discard_untrimmed: bool,
    ) -> Self {
        let mut bins = HashMap::new();

        // All (name1, name2) combinations.
        for n1 in adapter_names1 {
            for n2 in adapter_names2 {
                bins.insert((Some(n1.clone()), Some(n2.clone())), Vec::new());
            }
        }

        // If not discarding untrimmed, add extra bins for unmatched reads.
        if !discard_untrimmed {
            bins.insert((None, None), Vec::new());
            for n2 in adapter_names2 {
                bins.insert((None, Some(n2.clone())), Vec::new());
            }
            for n1 in adapter_names1 {
                bins.insert((Some(n1.clone()), None), Vec::new());
            }
        }

        Self {
            bins,
            statistics: ReadLengthStatistics::new(),
        }
    }
}

impl PairedEndStep for CombinatorialDemultiplexer {
    fn process_pair(
        &mut self,
        read1: SequenceRecord,
        read2: SequenceRecord,
        info1: &ModificationInfo,
        info2: &ModificationInfo,
    ) -> Option<RecordPair> {
        let name1 = if info1.matches.is_empty() {
            None
        } else {
            info1
                .matches
                .last()
                .and_then(|m| m.adapter_name.clone())
        };
        let name2 = if info2.matches.is_empty() {
            None
        } else {
            info2
                .matches
                .last()
                .and_then(|m| m.adapter_name.clone())
        };
        let key = (name1, name2);
        if let Some(bin) = self.bins.get_mut(&key) {
            self.statistics.update2(read1.len(), read2.len());
            bin.push((read1, read2));
        }
        None
    }
}

impl HasStatistics for CombinatorialDemultiplexer {
    fn get_statistics(&self) -> &ReadLengthStatistics {
        &self.statistics
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::info::{MatchInfo, ModificationInfo};
    use crate::predicates::TooShort;
    use crate::record::SequenceRecord;

    fn make_read(name: &str, seq: &str) -> SequenceRecord {
        SequenceRecord::new(name, seq, Some(&"I".repeat(seq.len())))
    }

    fn make_info(read: &SequenceRecord) -> ModificationInfo {
        ModificationInfo::new(read.clone())
    }

    fn make_info_with_adapter(read: &SequenceRecord, adapter_name: &str) -> ModificationInfo {
        let mut info = ModificationInfo::new(read.clone());
        info.matches.push(MatchInfo {
            adapter_name: Some(adapter_name.to_string()),
        });
        info
    }

    // -- SingleEndFilter --

    #[test]
    fn test_single_end_filter_passes() {
        let mut filter = SingleEndFilter::new(Box::new(TooShort::new(3)), false);
        let read = make_read("r1", "ACGTACGT");
        let info = make_info(&read);
        let result = filter.process(read.clone(), &info);
        assert!(result.is_some());
        assert_eq!(result.unwrap().sequence, "ACGTACGT");
        assert_eq!(filter.filtered(), 0);
    }

    #[test]
    fn test_single_end_filter_filters() {
        let mut filter = SingleEndFilter::new(Box::new(TooShort::new(10)), true);
        let read = make_read("r1", "ACGT");
        let info = make_info(&read);
        let result = filter.process(read, &info);
        assert!(result.is_none());
        assert_eq!(filter.filtered(), 1);

        let collected = filter.take_filtered_reads();
        assert_eq!(collected.len(), 1);
        assert_eq!(collected[0].sequence, "ACGT");
    }

    #[test]
    fn test_single_end_filter_identifier() {
        let filter = SingleEndFilter::new(Box::new(TooShort::new(10)), false);
        assert_eq!(filter.descriptive_identifier(), "too_short");
    }

    // -- PairedEndFilter --

    #[test]
    fn test_paired_end_filter_any_mode() {
        let mut filter = PairedEndFilter::new(
            Some(Box::new(TooShort::new(10))),
            Some(Box::new(TooShort::new(10))),
            PairFilterMode::Any,
            false,
        );
        let r1 = make_read("r1", "ACGT"); // too short
        let r2 = make_read("r2", "ACGTACGTACGT"); // long enough
        let info1 = make_info(&r1);
        let info2 = make_info(&r2);
        let result = filter.process_pair(r1, r2, &info1, &info2);
        assert!(result.is_none()); // filtered because r1 is too short
        assert_eq!(filter.filtered(), 1);
    }

    #[test]
    fn test_paired_end_filter_both_mode() {
        let mut filter = PairedEndFilter::new(
            Some(Box::new(TooShort::new(10))),
            Some(Box::new(TooShort::new(10))),
            PairFilterMode::Both,
            false,
        );
        let r1 = make_read("r1", "ACGT"); // too short
        let r2 = make_read("r2", "ACGTACGTACGT"); // long enough
        let info1 = make_info(&r1);
        let info2 = make_info(&r2);
        let result = filter.process_pair(r1, r2, &info1, &info2);
        assert!(result.is_some()); // NOT filtered: both must match
        assert_eq!(filter.filtered(), 0);
    }

    #[test]
    fn test_paired_end_filter_both_mode_both_short() {
        let mut filter = PairedEndFilter::new(
            Some(Box::new(TooShort::new(10))),
            Some(Box::new(TooShort::new(10))),
            PairFilterMode::Both,
            false,
        );
        let r1 = make_read("r1", "ACGT");
        let r2 = make_read("r2", "TGCA");
        let info1 = make_info(&r1);
        let info2 = make_info(&r2);
        let result = filter.process_pair(r1, r2, &info1, &info2);
        assert!(result.is_none());
        assert_eq!(filter.filtered(), 1);
    }

    #[test]
    fn test_paired_end_filter_first_mode() {
        let mut filter = PairedEndFilter::new(
            Some(Box::new(TooShort::new(10))),
            Some(Box::new(TooShort::new(10))),
            PairFilterMode::First,
            false,
        );
        let r1 = make_read("r1", "ACGTACGTACGT"); // long enough
        let r2 = make_read("r2", "ACGT"); // too short, but irrelevant
        let info1 = make_info(&r1);
        let info2 = make_info(&r2);
        let result = filter.process_pair(r1, r2, &info1, &info2);
        assert!(result.is_some()); // not filtered: first read is ok
    }

    #[test]
    fn test_paired_end_filter_predicate1_only() {
        let mut filter = PairedEndFilter::new(
            Some(Box::new(TooShort::new(10))),
            None,
            PairFilterMode::Any,
            false,
        );
        let r1 = make_read("r1", "ACGT");
        let r2 = make_read("r2", "ACGTACGTACGT");
        let info1 = make_info(&r1);
        let info2 = make_info(&r2);
        let result = filter.process_pair(r1, r2, &info1, &info2);
        assert!(result.is_none()); // filtered by predicate1 only
    }

    #[test]
    fn test_paired_end_filter_collect() {
        let mut filter = PairedEndFilter::new(
            Some(Box::new(TooShort::new(10))),
            Some(Box::new(TooShort::new(10))),
            PairFilterMode::Any,
            true,
        );
        let r1 = make_read("r1", "ACGT");
        let r2 = make_read("r2", "TG");
        let info1 = make_info(&r1);
        let info2 = make_info(&r2);
        filter.process_pair(r1, r2, &info1, &info2);

        let pairs = filter.take_filtered_pairs();
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].0.name, "r1");
        assert_eq!(pairs[0].1.name, "r2");
    }

    // -- SingleEndSink --

    #[test]
    fn test_single_end_sink() {
        let mut sink = SingleEndSink::new();
        let r1 = make_read("r1", "ACGT");
        let r2 = make_read("r2", "TGCATGCA");
        let info1 = make_info(&r1);
        let info2 = make_info(&r2);

        assert!(sink.process(r1, &info1).is_none());
        assert!(sink.process(r2, &info2).is_none());

        assert_eq!(sink.reads.len(), 2);
        assert_eq!(sink.get_statistics().written_reads(), 2);
        let (bp1, _bp2) = sink.get_statistics().written_bp();
        assert_eq!(bp1, 12); // 4 + 8
    }

    // -- PairedEndSink --

    #[test]
    fn test_paired_end_sink() {
        let mut sink = PairedEndSink::new();
        let r1 = make_read("r1", "ACGT");
        let r2 = make_read("r2", "TGCATG");
        let info1 = make_info(&r1);
        let info2 = make_info(&r2);

        assert!(sink.process_pair(r1, r2, &info1, &info2).is_none());

        assert_eq!(sink.pairs.len(), 1);
        assert_eq!(sink.get_statistics().written_reads(), 1);
        let (bp1, bp2) = sink.get_statistics().written_bp();
        assert_eq!(bp1, 4);
        assert_eq!(bp2, 6);
    }

    // -- Demultiplexer --

    #[test]
    fn test_demultiplexer_routes_by_adapter() {
        let names = vec!["adapter_a".to_string(), "adapter_b".to_string()];
        let mut demux = Demultiplexer::new(&names, false);

        let r1 = make_read("r1", "ACGT");
        let info1 = make_info_with_adapter(&r1, "adapter_a");
        demux.process(r1, &info1);

        let r2 = make_read("r2", "TGCA");
        let info2 = make_info_with_adapter(&r2, "adapter_b");
        demux.process(r2, &info2);

        assert_eq!(demux.bins["adapter_a"].len(), 1);
        assert_eq!(demux.bins["adapter_b"].len(), 1);
        assert_eq!(demux.bins["adapter_a"][0].name, "r1");
        assert_eq!(demux.bins["adapter_b"][0].name, "r2");
    }

    #[test]
    fn test_demultiplexer_untrimmed() {
        let names = vec!["adapter_a".to_string()];
        let mut demux = Demultiplexer::new(&names, false);

        let r1 = make_read("r1", "ACGT");
        let info1 = make_info(&r1); // no matches
        demux.process(r1, &info1);

        assert_eq!(demux.untrimmed.as_ref().unwrap().len(), 1);
        assert_eq!(demux.filtered(), 0);
    }

    #[test]
    fn test_demultiplexer_discard_untrimmed() {
        let names = vec!["adapter_a".to_string()];
        let mut demux = Demultiplexer::new(&names, true);

        let r1 = make_read("r1", "ACGT");
        let info1 = make_info(&r1);
        demux.process(r1, &info1);

        assert!(demux.untrimmed.is_none());
        assert_eq!(demux.filtered(), 1);
    }

    // -- PairedDemultiplexer --

    #[test]
    fn test_paired_demultiplexer() {
        let names = vec!["adapter_a".to_string()];
        let mut demux = PairedDemultiplexer::new(&names, false);

        let r1 = make_read("r1", "ACGT");
        let r2 = make_read("r2", "TGCA");
        let info1 = make_info_with_adapter(&r1, "adapter_a");
        let info2 = make_info(&r2);
        demux.process_pair(r1, r2, &info1, &info2);

        assert_eq!(demux.bins["adapter_a"].len(), 1);
    }

    #[test]
    fn test_paired_demultiplexer_untrimmed() {
        let names = vec!["adapter_a".to_string()];
        let mut demux = PairedDemultiplexer::new(&names, false);

        let r1 = make_read("r1", "ACGT");
        let r2 = make_read("r2", "TGCA");
        let info1 = make_info(&r1);
        let info2 = make_info(&r2);
        demux.process_pair(r1, r2, &info1, &info2);

        assert_eq!(demux.untrimmed.as_ref().unwrap().len(), 1);
    }

    // -- CombinatorialDemultiplexer --

    #[test]
    fn test_combinatorial_demux_routes() {
        let names1 = vec!["a1".to_string()];
        let names2 = vec!["b1".to_string()];
        let mut demux = CombinatorialDemultiplexer::new(&names1, &names2, false);

        let r1 = make_read("r1", "ACGT");
        let r2 = make_read("r2", "TGCA");
        let info1 = make_info_with_adapter(&r1, "a1");
        let info2 = make_info_with_adapter(&r2, "b1");
        demux.process_pair(r1, r2, &info1, &info2);

        let key = (Some("a1".to_string()), Some("b1".to_string()));
        assert_eq!(demux.bins[&key].len(), 1);
    }

    #[test]
    fn test_combinatorial_demux_untrimmed_bins() {
        let names1 = vec!["a1".to_string()];
        let names2 = vec!["b1".to_string()];
        let demux = CombinatorialDemultiplexer::new(&names1, &names2, false);

        // Should have bins for all combinations including None keys.
        assert!(demux.bins.contains_key(&(Some("a1".to_string()), Some("b1".to_string()))));
        assert!(demux.bins.contains_key(&(None, None)));
        assert!(demux.bins.contains_key(&(None, Some("b1".to_string()))));
        assert!(demux.bins.contains_key(&(Some("a1".to_string()), None)));
    }

    #[test]
    fn test_combinatorial_demux_discard_untrimmed() {
        let names1 = vec!["a1".to_string()];
        let names2 = vec!["b1".to_string()];
        let demux = CombinatorialDemultiplexer::new(&names1, &names2, true);

        // Should only have the fully-matched bin.
        assert!(demux.bins.contains_key(&(Some("a1".to_string()), Some("b1".to_string()))));
        assert!(!demux.bins.contains_key(&(None, None)));
    }

    #[test]
    fn test_combinatorial_demux_no_match_discarded() {
        let names1 = vec!["a1".to_string()];
        let names2 = vec!["b1".to_string()];
        let mut demux = CombinatorialDemultiplexer::new(&names1, &names2, true);

        let r1 = make_read("r1", "ACGT");
        let r2 = make_read("r2", "TGCA");
        let info1 = make_info(&r1); // no match
        let info2 = make_info(&r2); // no match
        demux.process_pair(r1, r2, &info1, &info2);

        // (None, None) bin does not exist when discarding untrimmed.
        assert_eq!(demux.get_statistics().written_reads(), 0);
    }
}
