// Adapter finding and trimming module.
//
// Port of cutadapt/adapters.py.
//
// This module provides adapter types for finding and trimming adapters
// in sequencing reads. The main types are:
//
// - `SingleAdapter` enum: wraps 9 concrete adapter types
// - `LinkedAdapter`: a front adapter combined with a back adapter
// - `MultipleAdapters`: tries all adapters, keeps the best match
// - `AdapterMatch` enum: the match result (RemoveBefore, RemoveAfter, Linked)

pub mod linked;
pub mod matches;
pub mod single;
pub mod statistics;

use std::sync::atomic::{AtomicU64, Ordering};

use crate::align::EndSkip;

pub use linked::LinkedAdapter;
pub use matches::{AdapterMatch, LinkedMatch, RemoveAfterMatch, RemoveBeforeMatch};
pub use single::{
    AnywhereAdapter, BackAdapter, FrontAdapter, NonInternalBackAdapter, NonInternalFrontAdapter,
    PrefixAdapter, RightmostBackAdapter, RightmostFrontAdapter, SingleAdapter,
    SingleAdapterParams, SuffixAdapter,
};
pub use statistics::{
    AdapterStatistics, AnywhereAdapterStatistics, BackAdapterStatistics, EndStatistics,
    FrontAdapterStatistics, LinkedAdapterStatistics,
};

// ---------------------------------------------------------------------------
// Where enum: aligner flag combinations for each adapter type
// ---------------------------------------------------------------------------

/// Aligner flag combinations for adapter types.
///
/// "REFERENCE" is the adapter sequence, "QUERY" is the read sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Where(u8);

impl Where {
    pub const BACK: Self = Self(
        EndSkip::QUERY_START.bits() | EndSkip::QUERY_STOP.bits() | EndSkip::REFERENCE_END.bits(),
    ); // 14

    pub const FRONT: Self = Self(
        EndSkip::QUERY_START.bits() | EndSkip::QUERY_STOP.bits() | EndSkip::REFERENCE_START.bits(),
    ); // 11

    pub const PREFIX: Self = Self(EndSkip::QUERY_STOP.bits()); // 8

    pub const SUFFIX: Self = Self(EndSkip::QUERY_START.bits()); // 2

    pub const FRONT_NOT_INTERNAL: Self =
        Self(EndSkip::REFERENCE_START.bits() | EndSkip::QUERY_STOP.bits()); // 9

    pub const BACK_NOT_INTERNAL: Self =
        Self(EndSkip::QUERY_START.bits() | EndSkip::REFERENCE_END.bits()); // 6

    pub const ANYWHERE: Self = Self(EndSkip::SEMIGLOBAL.bits()); // 15

    pub fn bits(self) -> u8 {
        self.0
    }
}

// ---------------------------------------------------------------------------
// Global adapter name counter
// ---------------------------------------------------------------------------

static ADAPTER_NAME_COUNTER: AtomicU64 = AtomicU64::new(1);

/// Generate a unique adapter name.
pub fn generate_adapter_name() -> String {
    ADAPTER_NAME_COUNTER
        .fetch_add(1, Ordering::Relaxed)
        .to_string()
}

/// Reset the adapter name counter (primarily for testing).
pub fn reset_adapter_name_counter() {
    ADAPTER_NAME_COUNTER.store(1, Ordering::Relaxed);
}

// ---------------------------------------------------------------------------
// Adapter enum (wrapping single and linked)
// ---------------------------------------------------------------------------

/// Top-level adapter enum.
#[derive(Debug)]
pub enum Adapter {
    Single(SingleAdapter),
    Linked(LinkedAdapter),
}

impl Adapter {
    pub fn name(&self) -> &str {
        match self {
            Adapter::Single(a) => a.name(),
            Adapter::Linked(a) => &a.name,
        }
    }

    pub fn match_to(&mut self, sequence: &str) -> Option<AdapterMatch> {
        match self {
            Adapter::Single(a) => a.match_to(sequence),
            Adapter::Linked(a) => a.match_to(sequence),
        }
    }

    pub fn create_statistics(&self) -> Box<dyn AdapterStatistics> {
        match self {
            Adapter::Single(a) => a.create_statistics(),
            Adapter::Linked(a) => a.create_statistics(),
        }
    }

    pub fn descriptive_identifier(&self) -> &str {
        match self {
            Adapter::Single(a) => a.descriptive_identifier(),
            Adapter::Linked(a) => a.descriptive_identifier(),
        }
    }

    pub fn spec(&self) -> String {
        match self {
            Adapter::Single(a) => a.spec(),
            Adapter::Linked(a) => a.spec(),
        }
    }
}

// ---------------------------------------------------------------------------
// MultipleAdapters
// ---------------------------------------------------------------------------

/// Tries multiple adapters and returns the best match.
#[derive(Debug)]
pub struct MultipleAdapters {
    pub adapters: Vec<Adapter>,
}

impl MultipleAdapters {
    pub fn new(adapters: Vec<Adapter>) -> Self {
        Self { adapters }
    }

    /// Find the adapter that best matches the sequence.
    /// Returns the match and the index of the adapter that produced it.
    pub fn match_to(&mut self, sequence: &str) -> Option<(AdapterMatch, usize)> {
        let mut best_match: Option<AdapterMatch> = None;
        let mut best_index = 0;

        for (i, adapter) in self.adapters.iter_mut().enumerate() {
            if let Some(m) = adapter.match_to(sequence) {
                let dominated = if let Some(ref bm) = best_match {
                    m.score() > bm.score()
                        || (m.score() == bm.score() && m.errors() < bm.errors())
                } else {
                    true
                };
                if dominated {
                    best_match = Some(m);
                    best_index = i;
                }
            }
        }
        best_match.map(|m| (m, best_index))
    }

    pub fn len(&self) -> usize {
        self.adapters.len()
    }

    pub fn is_empty(&self) -> bool {
        self.adapters.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    static ADAPTER_NAME_TEST_LOCK: Mutex<()> = Mutex::new(());

    #[test]
    fn test_where_values() {
        assert_eq!(Where::BACK.bits(), 14);
        assert_eq!(Where::FRONT.bits(), 11);
        assert_eq!(Where::PREFIX.bits(), 8);
        assert_eq!(Where::SUFFIX.bits(), 2);
        assert_eq!(Where::FRONT_NOT_INTERNAL.bits(), 9);
        assert_eq!(Where::BACK_NOT_INTERNAL.bits(), 6);
        assert_eq!(Where::ANYWHERE.bits(), 15);
    }

    #[test]
    fn test_generate_adapter_name() {
        let _guard = ADAPTER_NAME_TEST_LOCK.lock().unwrap();
        reset_adapter_name_counter();
        let name1 = generate_adapter_name();
        let name2 = generate_adapter_name();
        assert_eq!(name1, "1");
        assert_eq!(name2, "2");
    }

    #[test]
    fn test_adapter_enum_single() {
        let _guard = ADAPTER_NAME_TEST_LOCK.lock().unwrap();
        reset_adapter_name_counter();
        let params = SingleAdapterParams {
            sequence: "AGATCGGAAGAGC".to_string(),
            max_errors: 0.1,
            min_overlap: 3,
            read_wildcards: false,
            adapter_wildcards: true,
            name: Some("test".to_string()),
            indels: true,
        };
        let mut adapter = Adapter::Single(SingleAdapter::Back(BackAdapter::new(params)));
        assert_eq!(adapter.name(), "test");

        let m = adapter.match_to("ACGTACGTAGATCGGAAGAGC");
        assert!(m.is_some());
    }

    #[test]
    fn test_multiple_adapters() {
        let _guard = ADAPTER_NAME_TEST_LOCK.lock().unwrap();
        reset_adapter_name_counter();
        let params1 = SingleAdapterParams {
            sequence: "AAAA".to_string(),
            max_errors: 0.0,
            min_overlap: 4,
            read_wildcards: false,
            adapter_wildcards: false,
            name: Some("adapter_a".to_string()),
            indels: true,
        };
        let params2 = SingleAdapterParams {
            sequence: "CCCC".to_string(),
            max_errors: 0.0,
            min_overlap: 4,
            read_wildcards: false,
            adapter_wildcards: false,
            name: Some("adapter_c".to_string()),
            indels: true,
        };
        let adapters = vec![
            Adapter::Single(SingleAdapter::Back(BackAdapter::new(params1))),
            Adapter::Single(SingleAdapter::Back(BackAdapter::new(params2))),
        ];
        let mut multi = MultipleAdapters::new(adapters);
        assert_eq!(multi.len(), 2);

        let result = multi.match_to("ACGTACGTCCCC");
        assert!(result.is_some());
        let (_m, idx) = result.unwrap();
        assert_eq!(idx, 1); // CCCC adapter should match
    }

    #[test]
    fn test_multiple_adapters_no_match() {
        let _guard = ADAPTER_NAME_TEST_LOCK.lock().unwrap();
        reset_adapter_name_counter();
        let params = SingleAdapterParams {
            sequence: "AGATCGGAAGAGC".to_string(),
            max_errors: 0.1,
            min_overlap: 3,
            read_wildcards: false,
            adapter_wildcards: true,
            name: Some("test".to_string()),
            indels: true,
        };
        let adapters = vec![Adapter::Single(SingleAdapter::Back(BackAdapter::new(
            params,
        )))];
        let mut multi = MultipleAdapters::new(adapters);

        let m = multi.match_to("TTTTTTTTTTTTTTTTTTTTTT");
        assert!(m.is_none());
    }
}
