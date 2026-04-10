// Adapter statistics collection.
//
// Port of EndStatistics, AdapterStatistics, and related classes
// from cutadapt/adapters.py.

use std::collections::HashMap;

use super::matches::{AdapterMatch, LinkedMatch};

/// Statistics about matching at a single end (5' or 3').
#[derive(Debug, Clone)]
pub struct EndStatistics {
    pub max_error_rate: f64,
    pub sequence: String,
    pub effective_length: i32,
    pub has_wildcards: bool,
    pub indels: bool,
    pub adapter_type: String,
    pub allows_partial_matches: bool,
    /// errors[length][error_count] = number of times observed
    pub errors: HashMap<usize, HashMap<i32, u64>>,
    pub adjacent_bases: HashMap<String, u64>,
    /// Whether this end removes a prefix (true for front adapters).
    remove_prefix: bool,
}

impl EndStatistics {
    pub fn new(
        max_error_rate: f64,
        sequence: &str,
        effective_length: i32,
        has_wildcards: bool,
        indels: bool,
        adapter_type: &str,
        allows_partial_matches: bool,
        remove_prefix: bool,
    ) -> Self {
        let mut adjacent_bases = HashMap::new();
        for base in &["A", "C", "G", "T", ""] {
            adjacent_bases.insert(base.to_string(), 0);
        }
        Self {
            max_error_rate,
            sequence: sequence.to_string(),
            effective_length,
            has_wildcards,
            indels,
            adapter_type: adapter_type.to_string(),
            allows_partial_matches,
            errors: HashMap::new(),
            adjacent_bases,
            remove_prefix,
        }
    }

    /// Add error count for a given length.
    pub fn add_error(&mut self, length: usize, errors: i32) {
        *self
            .errors
            .entry(length)
            .or_default()
            .entry(errors)
            .or_insert(0) += 1;
    }

    /// Add an adjacent base observation.
    pub fn add_adjacent_base(&mut self, base: &str) {
        if let Some(count) = self.adjacent_bases.get_mut(base) {
            *count += 1;
        } else {
            *self.adjacent_bases.entry(String::new()).or_insert(0) += 1;
        }
    }

    /// Get a map from lengths to total number of matches at that length.
    pub fn lengths(&self) -> HashMap<usize, u64> {
        self.errors
            .iter()
            .map(|(&length, error_dict)| (length, error_dict.values().sum()))
            .collect()
    }

    /// Merge another EndStatistics into this one.
    pub fn merge(&mut self, other: &EndStatistics) {
        assert_eq!(self.max_error_rate, other.max_error_rate);
        assert_eq!(self.sequence, other.sequence);
        assert_eq!(self.effective_length, other.effective_length);
        assert_eq!(self.indels, other.indels);

        for base in &["A", "C", "G", "T", ""] {
            let key = base.to_string();
            *self.adjacent_bases.entry(key.clone()).or_insert(0) +=
                other.adjacent_bases.get(&key).copied().unwrap_or(0);
        }
        for (&length, error_dict) in &other.errors {
            for (&errors, &count) in error_dict {
                *self
                    .errors
                    .entry(length)
                    .or_default()
                    .entry(errors)
                    .or_insert(0) += count;
            }
        }
    }

    /// Estimate probabilities that this adapter end matches a random sequence.
    pub fn random_match_probabilities(&self, gc_content: f64) -> Vec<f64> {
        assert!((0.0..=1.0).contains(&gc_content));
        let seq: Vec<u8> = if self.remove_prefix {
            self.sequence.bytes().rev().collect()
        } else {
            self.sequence.bytes().collect()
        };
        let allowed_bytes: &[u8] = if self.has_wildcards {
            b"CGRYSKMBDHVN"
        } else {
            b"GC"
        };
        let mut p = 1.0f64;
        let mut probabilities = vec![p];
        for &c in &seq {
            if allowed_bytes.contains(&c) {
                p *= gc_content / 2.0;
            } else {
                p *= (1.0 - gc_content) / 2.0;
            }
            probabilities.push(p);
        }
        probabilities
    }
}

/// Statistics about a single adapter (abstract-ish interface).
pub trait AdapterStatistics: std::fmt::Debug + Send {
    fn name(&self) -> &str;
    fn reverse_complemented(&self) -> u64;
    fn add_reverse_complemented(&mut self, count: u64);
    fn add_match(&mut self, m: &AdapterMatch);
    fn end_statistics(&self) -> (Option<&EndStatistics>, Option<&EndStatistics>);
    fn merge(&mut self, other: &dyn AdapterStatistics);
}

/// Statistics for a front (5') adapter.
#[derive(Debug, Clone)]
pub struct FrontAdapterStatistics {
    pub name: String,
    pub end: EndStatistics,
    pub reverse_complemented_count: u64,
}

impl FrontAdapterStatistics {
    pub fn new(name: &str, end: EndStatistics) -> Self {
        Self {
            name: name.to_string(),
            end,
            reverse_complemented_count: 0,
        }
    }
}

impl AdapterStatistics for FrontAdapterStatistics {
    fn name(&self) -> &str {
        &self.name
    }
    fn reverse_complemented(&self) -> u64 {
        self.reverse_complemented_count
    }
    fn add_reverse_complemented(&mut self, count: u64) {
        self.reverse_complemented_count += count;
    }
    fn add_match(&mut self, m: &AdapterMatch) {
        if let AdapterMatch::RemoveBefore(ref rbm) = m {
            self.end
                .add_error(rbm.removed_sequence_length(), rbm.errors());
        }
    }
    fn end_statistics(&self) -> (Option<&EndStatistics>, Option<&EndStatistics>) {
        (Some(&self.end), None)
    }
    fn merge(&mut self, other: &dyn AdapterStatistics) {
        if let (Some(other_end), _) = other.end_statistics() {
            self.end.merge(other_end);
        }
        self.reverse_complemented_count += other.reverse_complemented();
    }
}

/// Statistics for a back (3') adapter.
#[derive(Debug, Clone)]
pub struct BackAdapterStatistics {
    pub name: String,
    pub end: EndStatistics,
    pub reverse_complemented_count: u64,
}

impl BackAdapterStatistics {
    pub fn new(name: &str, end: EndStatistics) -> Self {
        Self {
            name: name.to_string(),
            end,
            reverse_complemented_count: 0,
        }
    }
}

impl AdapterStatistics for BackAdapterStatistics {
    fn name(&self) -> &str {
        &self.name
    }
    fn reverse_complemented(&self) -> u64 {
        self.reverse_complemented_count
    }
    fn add_reverse_complemented(&mut self, count: u64) {
        self.reverse_complemented_count += count;
    }
    fn add_match(&mut self, m: &AdapterMatch) {
        if let AdapterMatch::RemoveAfter(ref ram) = m {
            self.end
                .add_error(ram.removed_sequence_length(), ram.errors());
            self.end.add_adjacent_base(ram.adjacent_base());
        }
    }
    fn end_statistics(&self) -> (Option<&EndStatistics>, Option<&EndStatistics>) {
        (None, Some(&self.end))
    }
    fn merge(&mut self, other: &dyn AdapterStatistics) {
        if let (_, Some(other_end)) = other.end_statistics() {
            self.end.merge(other_end);
        }
        self.reverse_complemented_count += other.reverse_complemented();
    }
}

/// Statistics for an anywhere adapter (can be 5' or 3').
#[derive(Debug, Clone)]
pub struct AnywhereAdapterStatistics {
    pub name: String,
    pub front: EndStatistics,
    pub back: EndStatistics,
    pub reverse_complemented_count: u64,
}

impl AnywhereAdapterStatistics {
    pub fn new(name: &str, front: EndStatistics, back: EndStatistics) -> Self {
        Self {
            name: name.to_string(),
            front,
            back,
            reverse_complemented_count: 0,
        }
    }
}

impl AdapterStatistics for AnywhereAdapterStatistics {
    fn name(&self) -> &str {
        &self.name
    }
    fn reverse_complemented(&self) -> u64 {
        self.reverse_complemented_count
    }
    fn add_reverse_complemented(&mut self, count: u64) {
        self.reverse_complemented_count += count;
    }
    fn add_match(&mut self, m: &AdapterMatch) {
        match m {
            AdapterMatch::RemoveBefore(ref rbm) => {
                self.front
                    .add_error(rbm.removed_sequence_length(), rbm.errors());
            }
            AdapterMatch::RemoveAfter(ref ram) => {
                self.back
                    .add_error(ram.removed_sequence_length(), ram.errors());
                self.back.add_adjacent_base(ram.adjacent_base());
            }
            _ => {}
        }
    }
    fn end_statistics(&self) -> (Option<&EndStatistics>, Option<&EndStatistics>) {
        (Some(&self.front), Some(&self.back))
    }
    fn merge(&mut self, other: &dyn AdapterStatistics) {
        let (other_front, other_back) = other.end_statistics();
        if let Some(of) = other_front {
            self.front.merge(of);
        }
        if let Some(ob) = other_back {
            self.back.merge(ob);
        }
        self.reverse_complemented_count += other.reverse_complemented();
    }
}

/// Statistics for a linked adapter.
#[derive(Debug, Clone)]
pub struct LinkedAdapterStatistics {
    pub name: String,
    pub front: EndStatistics,
    pub back: EndStatistics,
    pub reverse_complemented_count: u64,
}

impl LinkedAdapterStatistics {
    pub fn new(name: &str, front: EndStatistics, back: EndStatistics) -> Self {
        Self {
            name: name.to_string(),
            front,
            back,
            reverse_complemented_count: 0,
        }
    }
}

impl AdapterStatistics for LinkedAdapterStatistics {
    fn name(&self) -> &str {
        &self.name
    }
    fn reverse_complemented(&self) -> u64 {
        self.reverse_complemented_count
    }
    fn add_reverse_complemented(&mut self, count: u64) {
        self.reverse_complemented_count += count;
    }
    fn add_match(&mut self, m: &AdapterMatch) {
        if let AdapterMatch::Linked(ref lm) = m {
            self.add_linked_match(lm);
        }
    }
    fn end_statistics(&self) -> (Option<&EndStatistics>, Option<&EndStatistics>) {
        (Some(&self.front), Some(&self.back))
    }
    fn merge(&mut self, other: &dyn AdapterStatistics) {
        let (other_front, other_back) = other.end_statistics();
        if let Some(of) = other_front {
            self.front.merge(of);
        }
        if let Some(ob) = other_back {
            self.back.merge(ob);
        }
        self.reverse_complemented_count += other.reverse_complemented();
    }
}

impl LinkedAdapterStatistics {
    fn add_linked_match(&mut self, m: &LinkedMatch) {
        if let Some(ref fm) = m.front_match {
            self.front
                .add_error(fm.removed_sequence_length(), fm.errors());
        }
        if let Some(ref bm) = m.back_match {
            self.back
                .add_error(bm.removed_sequence_length(), bm.errors());
            self.back.add_adjacent_base(bm.adjacent_base());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_end_stats(remove_prefix: bool) -> EndStatistics {
        EndStatistics::new(0.1, "ACGTACGT", 8, false, true, "test", true, remove_prefix)
    }

    #[test]
    fn test_end_statistics_add_error() {
        let mut es = make_end_stats(false);
        es.add_error(5, 1);
        es.add_error(5, 1);
        es.add_error(5, 0);
        assert_eq!(*es.errors.get(&5).unwrap().get(&1).unwrap(), 2);
        assert_eq!(*es.errors.get(&5).unwrap().get(&0).unwrap(), 1);
    }

    #[test]
    fn test_end_statistics_lengths() {
        let mut es = make_end_stats(false);
        es.add_error(5, 1);
        es.add_error(5, 0);
        es.add_error(8, 2);
        let lengths = es.lengths();
        assert_eq!(*lengths.get(&5).unwrap(), 2);
        assert_eq!(*lengths.get(&8).unwrap(), 1);
    }

    #[test]
    fn test_end_statistics_adjacent_base() {
        let mut es = make_end_stats(false);
        es.add_adjacent_base("A");
        es.add_adjacent_base("A");
        es.add_adjacent_base("C");
        es.add_adjacent_base("X"); // unknown -> goes to ""
        assert_eq!(*es.adjacent_bases.get("A").unwrap(), 2);
        assert_eq!(*es.adjacent_bases.get("C").unwrap(), 1);
        assert_eq!(*es.adjacent_bases.get("").unwrap(), 1);
    }

    #[test]
    fn test_end_statistics_merge() {
        let mut es1 = make_end_stats(false);
        es1.add_error(5, 1);
        let mut es2 = make_end_stats(false);
        es2.add_error(5, 1);
        es2.add_error(8, 0);
        es1.merge(&es2);
        assert_eq!(*es1.errors.get(&5).unwrap().get(&1).unwrap(), 2);
        assert_eq!(*es1.errors.get(&8).unwrap().get(&0).unwrap(), 1);
    }

    #[test]
    fn test_random_match_probabilities() {
        let es = make_end_stats(false);
        let probs = es.random_match_probabilities(0.5);
        assert_eq!(probs.len(), 9); // sequence length + 1
        assert!((probs[0] - 1.0).abs() < 1e-10);
        assert!(probs[1] < 1.0);
    }

    #[test]
    fn test_front_adapter_statistics() {
        let es = make_end_stats(true);
        let mut stats = FrontAdapterStatistics::new("test_adapter", es);
        let m = AdapterMatch::RemoveBefore(super::super::matches::RemoveBeforeMatch::new(
            0, 5, 0, 5, 5, 1, "ACGTACGT", "test_adapter", "ACGTA",
        ));
        stats.add_match(&m);
        let (front, back) = stats.end_statistics();
        assert!(front.is_some());
        assert!(back.is_none());
    }

    #[test]
    fn test_back_adapter_statistics() {
        let es = make_end_stats(false);
        let mut stats = BackAdapterStatistics::new("test_adapter", es);
        let m = AdapterMatch::RemoveAfter(super::super::matches::RemoveAfterMatch::new(
            0, 5, 3, 8, 5, 1, "ACGTACGT", "test_adapter", "ACGTA",
        ));
        stats.add_match(&m);
        let (front, back) = stats.end_statistics();
        assert!(front.is_none());
        assert!(back.is_some());
    }
}
