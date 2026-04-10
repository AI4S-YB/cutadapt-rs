// Single adapter types.
//
// Port of the SingleAdapter class hierarchy from cutadapt/adapters.py.
// Each adapter type differs in: aligner flags, kmer_finder parameters,
// match return type, and for Rightmost variants the reference is reversed.

use crate::align::{AlignmentResult, Aligner, PrefixComparer, SuffixComparer};
use crate::kmer_finder::KmerFinder;
use crate::kmer_heuristic::create_positions_and_kmers;

use super::matches::{AdapterMatch, RemoveAfterMatch, RemoveBeforeMatch};
use super::statistics::{
    AdapterStatistics, AnywhereAdapterStatistics, BackAdapterStatistics, EndStatistics,
    FrontAdapterStatistics,
};
use super::Where;

// ---------------------------------------------------------------------------
// Aligner wrapper enum (Aligner needs &mut self for locate, comparers need &self)
// ---------------------------------------------------------------------------

enum AlignerEnum {
    Full(Aligner),
    Prefix(PrefixComparer),
    Suffix(SuffixComparer),
}

impl AlignerEnum {
    fn locate(&mut self, query: &str) -> Option<AlignmentResult> {
        match self {
            AlignerEnum::Full(a) => a.locate(query),
            AlignerEnum::Prefix(c) => c.locate(query),
            AlignerEnum::Suffix(c) => c.locate(query),
        }
    }

    fn effective_length(&self) -> i32 {
        match self {
            AlignerEnum::Full(a) => a.effective_length,
            AlignerEnum::Prefix(c) => c.effective_length,
            AlignerEnum::Suffix(c) => c.effective_length(),
        }
    }
}

// ---------------------------------------------------------------------------
// Kmer finder wrapper (real vs. mock that always returns true)
// ---------------------------------------------------------------------------

enum KmerFinderEnum {
    Real(KmerFinder),
    Mock,
}

impl KmerFinderEnum {
    fn kmers_present(&self, sequence: &str) -> bool {
        match self {
            KmerFinderEnum::Real(kf) => kf.kmers_present(sequence),
            KmerFinderEnum::Mock => true,
        }
    }
}

// ---------------------------------------------------------------------------
// Common SingleAdapter fields
// ---------------------------------------------------------------------------

/// Configuration and state common to all single adapter types.
#[allow(dead_code)]
struct SingleAdapterCore {
    name: String,
    sequence: String,
    max_error_rate: f64,
    min_overlap: usize,
    adapter_wildcards: bool,
    read_wildcards: bool,
    indels: bool,
    allows_partial_matches: bool,
    aligner: AlignerEnum,
    kmer_finder: KmerFinderEnum,
}

impl std::fmt::Debug for SingleAdapterCore {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SingleAdapterCore")
            .field("name", &self.name)
            .field("sequence", &self.sequence)
            .field("max_error_rate", &self.max_error_rate)
            .field("min_overlap", &self.min_overlap)
            .finish()
    }
}

/// Builder-style parameters for creating a SingleAdapter.
#[derive(Debug, Clone)]
pub struct SingleAdapterParams {
    pub sequence: String,
    pub max_errors: f64,
    pub min_overlap: usize,
    pub read_wildcards: bool,
    pub adapter_wildcards: bool,
    pub name: Option<String>,
    pub indels: bool,
}

impl Default for SingleAdapterParams {
    fn default() -> Self {
        Self {
            sequence: String::new(),
            max_errors: 0.1,
            min_overlap: 3,
            read_wildcards: false,
            adapter_wildcards: true,
            name: None,
            indels: true,
        }
    }
}

// Helper to normalise parameters like the Python __init__.
fn normalise_params(params: &SingleAdapterParams) -> (String, f64, usize, bool) {
    let sequence = params
        .sequence
        .to_uppercase()
        .replace('U', "T")
        .replace('I', "N");
    assert!(!sequence.is_empty(), "Adapter sequence is empty");

    let mut max_error_rate = params.max_errors;
    let non_n_count = sequence.len() - sequence.matches('N').count();
    if max_error_rate >= 1.0 && non_n_count > 0 {
        max_error_rate /= non_n_count as f64;
    }

    let min_overlap = params.min_overlap.min(sequence.len());

    // Only enable adapter wildcards if non-ACGT characters are present.
    let adapter_wildcards =
        params.adapter_wildcards && !sequence.chars().all(|c| "ACGT".contains(c));

    (sequence, max_error_rate, min_overlap, adapter_wildcards)
}

fn make_aligner(
    sequence: &str,
    flags: u8,
    max_error_rate: f64,
    adapter_wildcards: bool,
    read_wildcards: bool,
    indels: bool,
    min_overlap: usize,
) -> AlignerEnum {
    let indel_cost = if indels { 1 } else { 100_000 };
    AlignerEnum::Full(Aligner::new(
        sequence,
        max_error_rate,
        flags,
        adapter_wildcards,
        read_wildcards,
        indel_cost,
        min_overlap as i32,
    ))
}

fn make_kmer_finder(
    sequence: &str,
    min_overlap: usize,
    max_error_rate: f64,
    back_adapter: bool,
    front_adapter: bool,
    internal: bool,
    adapter_wildcards: bool,
    read_wildcards: bool,
) -> KmerFinderEnum {
    let positions_and_kmers = create_positions_and_kmers(
        sequence,
        min_overlap,
        max_error_rate,
        back_adapter,
        front_adapter,
        internal,
    );
    match KmerFinder::new(&positions_and_kmers, adapter_wildcards, read_wildcards) {
        kf => KmerFinderEnum::Real(kf),
    }
}

// ---------------------------------------------------------------------------
// The 10 adapter types as an enum for dispatch
// ---------------------------------------------------------------------------

/// All single adapter types.
#[derive(Debug)]
pub enum SingleAdapter {
    Front(FrontAdapter),
    Back(BackAdapter),
    Anywhere(AnywhereAdapter),
    Prefix(PrefixAdapter),
    Suffix(SuffixAdapter),
    NonInternalFront(NonInternalFrontAdapter),
    NonInternalBack(NonInternalBackAdapter),
    RightmostFront(RightmostFrontAdapter),
    RightmostBack(RightmostBackAdapter),
}

impl SingleAdapter {
    pub fn name(&self) -> &str {
        self.core().name.as_str()
    }

    pub fn sequence(&self) -> &str {
        self.core().sequence.as_str()
    }

    pub fn description(&self) -> &str {
        match self {
            SingleAdapter::Front(_) => "regular 5'",
            SingleAdapter::Back(_) => "regular 3'",
            SingleAdapter::Anywhere(_) => "variable 5'/3'",
            SingleAdapter::Prefix(_) => "anchored 5'",
            SingleAdapter::Suffix(_) => "anchored 3'",
            SingleAdapter::NonInternalFront(_) => "non-internal 5'",
            SingleAdapter::NonInternalBack(_) => "non-internal 3'",
            SingleAdapter::RightmostFront(_) => "rightmost 5'",
            SingleAdapter::RightmostBack(_) => "rightmost 3'",
        }
    }

    pub fn descriptive_identifier(&self) -> &str {
        match self {
            SingleAdapter::Front(_) => "regular_five_prime",
            SingleAdapter::Back(_) => "regular_three_prime",
            SingleAdapter::Anywhere(_) => "anywhere",
            SingleAdapter::Prefix(_) => "anchored_five_prime",
            SingleAdapter::Suffix(_) => "anchored_three_prime",
            SingleAdapter::NonInternalFront(_) => "noninternal_five_prime",
            SingleAdapter::NonInternalBack(_) => "noninternal_three_prime",
            SingleAdapter::RightmostFront(_) => "rightmost_five_prime",
            SingleAdapter::RightmostBack(_) => "rightmost_three_prime",
        }
    }

    pub fn spec(&self) -> String {
        let seq = self.core().sequence.as_str();
        match self {
            SingleAdapter::Front(_) => format!("{}...", seq),
            SingleAdapter::Back(_) => seq.to_string(),
            SingleAdapter::Anywhere(_) => format!("...{}...", seq),
            SingleAdapter::Prefix(_) => format!("^{}...", seq),
            SingleAdapter::Suffix(_) => format!("{}$", seq),
            SingleAdapter::NonInternalFront(_) => format!("X{}...", seq),
            SingleAdapter::NonInternalBack(_) => format!("{}X", seq),
            SingleAdapter::RightmostFront(_) => format!("{}...;rightmost", seq),
            SingleAdapter::RightmostBack(_) => format!("{};rightmost", seq),
        }
    }

    pub fn effective_length(&self) -> i32 {
        self.core().aligner.effective_length()
    }

    pub fn max_error_rate(&self) -> f64 {
        self.core().max_error_rate
    }

    pub fn adapter_wildcards(&self) -> bool {
        self.core().adapter_wildcards
    }

    pub fn indels(&self) -> bool {
        self.core().indels
    }

    pub fn allows_partial_matches(&self) -> bool {
        self.core().allows_partial_matches
    }

    pub fn match_to(&mut self, sequence: &str) -> Option<AdapterMatch> {
        match self {
            SingleAdapter::Front(a) => a.match_to(sequence),
            SingleAdapter::Back(a) => a.match_to(sequence),
            SingleAdapter::Anywhere(a) => a.match_to(sequence),
            SingleAdapter::Prefix(a) => a.match_to(sequence),
            SingleAdapter::Suffix(a) => a.match_to(sequence),
            SingleAdapter::NonInternalFront(a) => a.match_to(sequence),
            SingleAdapter::NonInternalBack(a) => a.match_to(sequence),
            SingleAdapter::RightmostFront(a) => a.match_to(sequence),
            SingleAdapter::RightmostBack(a) => a.match_to(sequence),
        }
    }

    pub fn create_statistics(&self) -> Box<dyn AdapterStatistics> {
        match self {
            SingleAdapter::Front(_)
            | SingleAdapter::NonInternalFront(_)
            | SingleAdapter::Prefix(_)
            | SingleAdapter::RightmostFront(_) => {
                let es = self.make_end_stats(true);
                Box::new(FrontAdapterStatistics::new(self.name(), es))
            }
            SingleAdapter::Back(_)
            | SingleAdapter::NonInternalBack(_)
            | SingleAdapter::Suffix(_)
            | SingleAdapter::RightmostBack(_) => {
                let es = self.make_end_stats(false);
                Box::new(BackAdapterStatistics::new(self.name(), es))
            }
            SingleAdapter::Anywhere(_) => {
                let front_es = self.make_end_stats(true);
                let back_es = self.make_end_stats(false);
                Box::new(AnywhereAdapterStatistics::new(
                    self.name(),
                    front_es,
                    back_es,
                ))
            }
        }
    }

    fn core(&self) -> &SingleAdapterCore {
        match self {
            SingleAdapter::Front(a) => &a.core,
            SingleAdapter::Back(a) => &a.core,
            SingleAdapter::Anywhere(a) => &a.core,
            SingleAdapter::Prefix(a) => &a.core,
            SingleAdapter::Suffix(a) => &a.core,
            SingleAdapter::NonInternalFront(a) => &a.core,
            SingleAdapter::NonInternalBack(a) => &a.core,
            SingleAdapter::RightmostFront(a) => &a.core,
            SingleAdapter::RightmostBack(a) => &a.core,
        }
    }

    fn make_end_stats(&self, remove_prefix: bool) -> EndStatistics {
        let c = self.core();
        EndStatistics::new(
            c.max_error_rate,
            &c.sequence,
            c.aligner.effective_length(),
            c.adapter_wildcards,
            c.indels,
            self.descriptive_identifier(),
            c.allows_partial_matches,
            remove_prefix,
        )
    }
}

// ---------------------------------------------------------------------------
// Helper: build alignment result as a match
// ---------------------------------------------------------------------------

fn alignment_to_remove_before(
    alignment: AlignmentResult,
    sequence: &str,
    adapter_name: &str,
    adapter_sequence: &str,
) -> AdapterMatch {
    let (astart, astop, rstart, rstop, score, errors) = alignment;
    AdapterMatch::RemoveBefore(RemoveBeforeMatch::new(
        astart,
        astop,
        rstart,
        rstop,
        score,
        errors,
        sequence,
        adapter_name,
        adapter_sequence,
    ))
}

fn alignment_to_remove_after(
    alignment: AlignmentResult,
    sequence: &str,
    adapter_name: &str,
    adapter_sequence: &str,
) -> AdapterMatch {
    let (astart, astop, rstart, rstop, score, errors) = alignment;
    AdapterMatch::RemoveAfter(RemoveAfterMatch::new(
        astart,
        astop,
        rstart,
        rstop,
        score,
        errors,
        sequence,
        adapter_name,
        adapter_sequence,
    ))
}

// ---------------------------------------------------------------------------
// FrontAdapter
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct FrontAdapter {
    core: SingleAdapterCore,
    #[allow(dead_code)]
    force_anywhere: bool,
}

impl FrontAdapter {
    pub fn new(params: SingleAdapterParams) -> Self {
        Self::with_force_anywhere(params, false)
    }

    pub fn with_force_anywhere(params: SingleAdapterParams, force_anywhere: bool) -> Self {
        let (sequence, max_error_rate, min_overlap, adapter_wildcards) = normalise_params(&params);
        let name = params
            .name
            .unwrap_or_else(|| super::generate_adapter_name());

        let flags = if force_anywhere {
            Where::ANYWHERE
        } else {
            Where::FRONT
        };
        let aligner = make_aligner(
            &sequence,
            flags.bits(),
            max_error_rate,
            adapter_wildcards,
            params.read_wildcards,
            params.indels,
            min_overlap,
        );
        let kmer_finder = make_kmer_finder(
            &sequence,
            min_overlap,
            max_error_rate,
            force_anywhere, // back_adapter
            true,           // front_adapter
            true,           // internal
            adapter_wildcards,
            params.read_wildcards,
        );

        Self {
            core: SingleAdapterCore {
                name,
                sequence,
                max_error_rate,
                min_overlap,
                adapter_wildcards,
                read_wildcards: params.read_wildcards,
                indels: params.indels,
                allows_partial_matches: true,
                aligner,
                kmer_finder,
            },
            force_anywhere,
        }
    }

    fn match_to(&mut self, sequence: &str) -> Option<AdapterMatch> {
        if !self.core.kmer_finder.kmers_present(sequence) {
            return None;
        }
        let alignment = self.core.aligner.locate(sequence)?;
        Some(alignment_to_remove_before(
            alignment,
            sequence,
            &self.core.name,
            &self.core.sequence,
        ))
    }
}

// ---------------------------------------------------------------------------
// RightmostFrontAdapter
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct RightmostFrontAdapter {
    core: SingleAdapterCore,
    #[allow(dead_code)]
    force_anywhere: bool,
    /// Original (non-reversed) sequence, for coordinate flipping.
    original_sequence: String,
}

impl RightmostFrontAdapter {
    pub fn new(params: SingleAdapterParams) -> Self {
        Self::with_force_anywhere(params, false)
    }

    pub fn with_force_anywhere(params: SingleAdapterParams, force_anywhere: bool) -> Self {
        let (sequence, max_error_rate, min_overlap, adapter_wildcards) = normalise_params(&params);
        let name = params
            .name
            .unwrap_or_else(|| super::generate_adapter_name());

        let reversed: String = sequence.chars().rev().collect();
        let flags = if force_anywhere {
            Where::ANYWHERE
        } else {
            Where::BACK
        };
        let aligner = make_aligner(
            &reversed,
            flags.bits(),
            max_error_rate,
            adapter_wildcards,
            params.read_wildcards,
            params.indels,
            min_overlap,
        );
        let kmer_finder = make_kmer_finder(
            &reversed,
            min_overlap,
            max_error_rate,
            true,           // back_adapter
            force_anywhere, // front_adapter
            true,           // internal
            adapter_wildcards,
            params.read_wildcards,
        );

        Self {
            original_sequence: sequence.clone(),
            core: SingleAdapterCore {
                name,
                sequence,
                max_error_rate,
                min_overlap,
                adapter_wildcards,
                read_wildcards: params.read_wildcards,
                indels: params.indels,
                allows_partial_matches: true,
                aligner,
                kmer_finder,
            },
            force_anywhere,
        }
    }

    fn match_to(&mut self, sequence: &str) -> Option<AdapterMatch> {
        let reversed_sequence: String = sequence.chars().rev().collect();
        if !self.core.kmer_finder.kmers_present(&reversed_sequence) {
            return None;
        }
        let alignment = self.core.aligner.locate(&reversed_sequence)?;
        let (ref_start, ref_end, query_start, query_end, score, errors) = alignment;
        let seq_len = self.original_sequence.len();
        let q_len = sequence.len();
        let flipped = (
            seq_len - ref_end,
            seq_len - ref_start,
            q_len - query_end,
            q_len - query_start,
            score,
            errors,
        );
        Some(alignment_to_remove_before(
            flipped,
            sequence,
            &self.core.name,
            &self.core.sequence,
        ))
    }
}

// ---------------------------------------------------------------------------
// BackAdapter
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct BackAdapter {
    core: SingleAdapterCore,
    #[allow(dead_code)]
    force_anywhere: bool,
}

impl BackAdapter {
    pub fn new(params: SingleAdapterParams) -> Self {
        Self::with_force_anywhere(params, false)
    }

    pub fn with_force_anywhere(params: SingleAdapterParams, force_anywhere: bool) -> Self {
        let (sequence, max_error_rate, min_overlap, adapter_wildcards) = normalise_params(&params);
        let name = params
            .name
            .unwrap_or_else(|| super::generate_adapter_name());

        let flags = if force_anywhere {
            Where::ANYWHERE
        } else {
            Where::BACK
        };
        let aligner = make_aligner(
            &sequence,
            flags.bits(),
            max_error_rate,
            adapter_wildcards,
            params.read_wildcards,
            params.indels,
            min_overlap,
        );
        let kmer_finder = make_kmer_finder(
            &sequence,
            min_overlap,
            max_error_rate,
            true,           // back_adapter
            force_anywhere, // front_adapter
            true,           // internal
            adapter_wildcards,
            params.read_wildcards,
        );

        Self {
            core: SingleAdapterCore {
                name,
                sequence,
                max_error_rate,
                min_overlap,
                adapter_wildcards,
                read_wildcards: params.read_wildcards,
                indels: params.indels,
                allows_partial_matches: true,
                aligner,
                kmer_finder,
            },
            force_anywhere,
        }
    }

    fn match_to(&mut self, sequence: &str) -> Option<AdapterMatch> {
        if !self.core.kmer_finder.kmers_present(sequence) {
            return None;
        }
        let alignment = self.core.aligner.locate(sequence)?;
        Some(alignment_to_remove_after(
            alignment,
            sequence,
            &self.core.name,
            &self.core.sequence,
        ))
    }
}

// ---------------------------------------------------------------------------
// RightmostBackAdapter
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct RightmostBackAdapter {
    core: SingleAdapterCore,
    #[allow(dead_code)]
    force_anywhere: bool,
    original_sequence: String,
}

impl RightmostBackAdapter {
    pub fn new(params: SingleAdapterParams) -> Self {
        Self::with_force_anywhere(params, false)
    }

    pub fn with_force_anywhere(params: SingleAdapterParams, force_anywhere: bool) -> Self {
        let (sequence, max_error_rate, min_overlap, adapter_wildcards) = normalise_params(&params);
        let name = params
            .name
            .unwrap_or_else(|| super::generate_adapter_name());

        let reversed: String = sequence.chars().rev().collect();
        let flags = if force_anywhere {
            Where::ANYWHERE
        } else {
            Where::FRONT
        };
        let aligner = make_aligner(
            &reversed,
            flags.bits(),
            max_error_rate,
            adapter_wildcards,
            params.read_wildcards,
            params.indels,
            min_overlap,
        );
        let kmer_finder = make_kmer_finder(
            &reversed,
            min_overlap,
            max_error_rate,
            force_anywhere, // back_adapter
            true,           // front_adapter
            true,           // internal
            adapter_wildcards,
            params.read_wildcards,
        );

        Self {
            original_sequence: sequence.clone(),
            core: SingleAdapterCore {
                name,
                sequence,
                max_error_rate,
                min_overlap,
                adapter_wildcards,
                read_wildcards: params.read_wildcards,
                indels: params.indels,
                allows_partial_matches: true,
                aligner,
                kmer_finder,
            },
            force_anywhere,
        }
    }

    fn match_to(&mut self, sequence: &str) -> Option<AdapterMatch> {
        let reversed_sequence: String = sequence.chars().rev().collect();
        if !self.core.kmer_finder.kmers_present(&reversed_sequence) {
            return None;
        }
        let alignment = self.core.aligner.locate(&reversed_sequence)?;
        let (ref_start, ref_end, query_start, query_end, score, errors) = alignment;
        let seq_len = self.original_sequence.len();
        let q_len = sequence.len();
        let flipped = (
            seq_len - ref_end,
            seq_len - ref_start,
            q_len - query_end,
            q_len - query_start,
            score,
            errors,
        );
        Some(alignment_to_remove_after(
            flipped,
            sequence,
            &self.core.name,
            &self.core.sequence,
        ))
    }
}

// ---------------------------------------------------------------------------
// AnywhereAdapter
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct AnywhereAdapter {
    core: SingleAdapterCore,
}

impl AnywhereAdapter {
    pub fn new(params: SingleAdapterParams) -> Self {
        let (sequence, max_error_rate, min_overlap, adapter_wildcards) = normalise_params(&params);
        let name = params
            .name
            .unwrap_or_else(|| super::generate_adapter_name());

        let aligner = make_aligner(
            &sequence,
            Where::ANYWHERE.bits(),
            max_error_rate,
            adapter_wildcards,
            params.read_wildcards,
            params.indels,
            min_overlap,
        );
        let kmer_finder = make_kmer_finder(
            &sequence,
            min_overlap,
            max_error_rate,
            true, // back_adapter
            true, // front_adapter
            true, // internal
            adapter_wildcards,
            params.read_wildcards,
        );

        Self {
            core: SingleAdapterCore {
                name,
                sequence,
                max_error_rate,
                min_overlap,
                adapter_wildcards,
                read_wildcards: params.read_wildcards,
                indels: params.indels,
                allows_partial_matches: true,
                aligner,
                kmer_finder,
            },
        }
    }

    fn match_to(&mut self, sequence: &str) -> Option<AdapterMatch> {
        if !self.core.kmer_finder.kmers_present(sequence) {
            return None;
        }
        let alignment = self.core.aligner.locate(sequence)?;
        // If alignment starts at pos 0 in the query, it's a 5' adapter.
        if alignment.2 == 0 {
            Some(alignment_to_remove_before(
                alignment,
                sequence,
                &self.core.name,
                &self.core.sequence,
            ))
        } else {
            Some(alignment_to_remove_after(
                alignment,
                sequence,
                &self.core.name,
                &self.core.sequence,
            ))
        }
    }
}

// ---------------------------------------------------------------------------
// NonInternalFrontAdapter
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct NonInternalFrontAdapter {
    core: SingleAdapterCore,
    #[allow(dead_code)]
    force_anywhere: bool,
}

impl NonInternalFrontAdapter {
    pub fn new(params: SingleAdapterParams) -> Self {
        Self::with_force_anywhere(params, false)
    }

    pub fn with_force_anywhere(params: SingleAdapterParams, force_anywhere: bool) -> Self {
        let (sequence, max_error_rate, min_overlap, adapter_wildcards) = normalise_params(&params);
        let name = params
            .name
            .unwrap_or_else(|| super::generate_adapter_name());

        let aligner = make_aligner(
            &sequence,
            Where::FRONT_NOT_INTERNAL.bits(),
            max_error_rate,
            adapter_wildcards,
            params.read_wildcards,
            params.indels,
            min_overlap,
        );
        let kmer_finder = make_kmer_finder(
            &sequence,
            min_overlap,
            max_error_rate,
            force_anywhere, // back_adapter
            true,           // front_adapter
            false,          // internal = false for non-internal
            adapter_wildcards,
            params.read_wildcards,
        );

        Self {
            core: SingleAdapterCore {
                name,
                sequence,
                max_error_rate,
                min_overlap,
                adapter_wildcards,
                read_wildcards: params.read_wildcards,
                indels: params.indels,
                allows_partial_matches: true,
                aligner,
                kmer_finder,
            },
            force_anywhere,
        }
    }

    fn match_to(&mut self, sequence: &str) -> Option<AdapterMatch> {
        if !self.core.kmer_finder.kmers_present(sequence) {
            return None;
        }
        let alignment = self.core.aligner.locate(sequence)?;
        Some(alignment_to_remove_before(
            alignment,
            sequence,
            &self.core.name,
            &self.core.sequence,
        ))
    }
}

// ---------------------------------------------------------------------------
// NonInternalBackAdapter
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct NonInternalBackAdapter {
    core: SingleAdapterCore,
    #[allow(dead_code)]
    force_anywhere: bool,
}

impl NonInternalBackAdapter {
    pub fn new(params: SingleAdapterParams) -> Self {
        Self::with_force_anywhere(params, false)
    }

    pub fn with_force_anywhere(params: SingleAdapterParams, force_anywhere: bool) -> Self {
        let (sequence, max_error_rate, min_overlap, adapter_wildcards) = normalise_params(&params);
        let name = params
            .name
            .unwrap_or_else(|| super::generate_adapter_name());

        let aligner = make_aligner(
            &sequence,
            Where::BACK_NOT_INTERNAL.bits(),
            max_error_rate,
            adapter_wildcards,
            params.read_wildcards,
            params.indels,
            min_overlap,
        );
        let kmer_finder = make_kmer_finder(
            &sequence,
            min_overlap,
            max_error_rate,
            true,           // back_adapter
            force_anywhere, // front_adapter
            false,          // internal = false for non-internal
            adapter_wildcards,
            params.read_wildcards,
        );

        Self {
            core: SingleAdapterCore {
                name,
                sequence,
                max_error_rate,
                min_overlap,
                adapter_wildcards,
                read_wildcards: params.read_wildcards,
                indels: params.indels,
                allows_partial_matches: true,
                aligner,
                kmer_finder,
            },
            force_anywhere,
        }
    }

    fn match_to(&mut self, sequence: &str) -> Option<AdapterMatch> {
        if !self.core.kmer_finder.kmers_present(sequence) {
            return None;
        }
        let alignment = self.core.aligner.locate(sequence)?;
        Some(alignment_to_remove_after(
            alignment,
            sequence,
            &self.core.name,
            &self.core.sequence,
        ))
    }
}

// ---------------------------------------------------------------------------
// PrefixAdapter (anchored 5')
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct PrefixAdapter {
    core: SingleAdapterCore,
    #[allow(dead_code)]
    force_anywhere: bool,
}

impl PrefixAdapter {
    pub fn new(params: SingleAdapterParams) -> Self {
        Self::with_force_anywhere(params, false)
    }

    pub fn with_force_anywhere(mut params: SingleAdapterParams, force_anywhere: bool) -> Self {
        // PrefixAdapter forces min_overlap = len(sequence)
        let seq_upper = params
            .sequence
            .to_uppercase()
            .replace('U', "T")
            .replace('I', "N");
        params.min_overlap = seq_upper.len();

        let (sequence, max_error_rate, min_overlap, adapter_wildcards) = normalise_params(&params);
        let name = params
            .name
            .unwrap_or_else(|| super::generate_adapter_name());

        let (aligner, use_prefix_comparer) = if !params.indels {
            (
                AlignerEnum::Prefix(PrefixComparer::new(
                    &sequence,
                    max_error_rate,
                    adapter_wildcards,
                    params.read_wildcards,
                    min_overlap,
                )),
                true,
            )
        } else {
            (
                make_aligner(
                    &sequence,
                    Where::PREFIX.bits(),
                    max_error_rate,
                    adapter_wildcards,
                    params.read_wildcards,
                    params.indels,
                    min_overlap,
                ),
                false,
            )
        };

        let kmer_finder = if use_prefix_comparer {
            KmerFinderEnum::Mock
        } else {
            make_kmer_finder(
                &sequence,
                min_overlap,
                max_error_rate,
                force_anywhere, // back_adapter
                true,           // front_adapter
                false,          // internal
                adapter_wildcards,
                params.read_wildcards,
            )
        };

        Self {
            core: SingleAdapterCore {
                name,
                sequence,
                max_error_rate,
                min_overlap,
                adapter_wildcards,
                read_wildcards: params.read_wildcards,
                indels: params.indels,
                allows_partial_matches: false,
                aligner,
                kmer_finder,
            },
            force_anywhere,
        }
    }

    fn match_to(&mut self, sequence: &str) -> Option<AdapterMatch> {
        if !self.core.kmer_finder.kmers_present(sequence) {
            return None;
        }
        let alignment = self.core.aligner.locate(sequence)?;
        Some(alignment_to_remove_before(
            alignment,
            sequence,
            &self.core.name,
            &self.core.sequence,
        ))
    }
}

// ---------------------------------------------------------------------------
// SuffixAdapter (anchored 3')
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct SuffixAdapter {
    core: SingleAdapterCore,
    #[allow(dead_code)]
    force_anywhere: bool,
}

impl SuffixAdapter {
    pub fn new(params: SingleAdapterParams) -> Self {
        Self::with_force_anywhere(params, false)
    }

    pub fn with_force_anywhere(mut params: SingleAdapterParams, force_anywhere: bool) -> Self {
        let seq_upper = params
            .sequence
            .to_uppercase()
            .replace('U', "T")
            .replace('I', "N");
        params.min_overlap = seq_upper.len();

        let (sequence, max_error_rate, min_overlap, adapter_wildcards) = normalise_params(&params);
        let name = params
            .name
            .unwrap_or_else(|| super::generate_adapter_name());

        let (aligner, use_suffix_comparer) = if !params.indels {
            (
                AlignerEnum::Suffix(SuffixComparer::new(
                    &sequence,
                    max_error_rate,
                    adapter_wildcards,
                    params.read_wildcards,
                    min_overlap,
                )),
                true,
            )
        } else {
            (
                make_aligner(
                    &sequence,
                    Where::SUFFIX.bits(),
                    max_error_rate,
                    adapter_wildcards,
                    params.read_wildcards,
                    params.indels,
                    min_overlap,
                ),
                false,
            )
        };

        let kmer_finder = if use_suffix_comparer {
            KmerFinderEnum::Mock
        } else {
            make_kmer_finder(
                &sequence,
                min_overlap,
                max_error_rate,
                true,           // back_adapter
                force_anywhere, // front_adapter
                false,          // internal
                adapter_wildcards,
                params.read_wildcards,
            )
        };

        Self {
            core: SingleAdapterCore {
                name,
                sequence,
                max_error_rate,
                min_overlap,
                adapter_wildcards,
                read_wildcards: params.read_wildcards,
                indels: params.indels,
                allows_partial_matches: false,
                aligner,
                kmer_finder,
            },
            force_anywhere,
        }
    }

    fn match_to(&mut self, sequence: &str) -> Option<AdapterMatch> {
        if !self.core.kmer_finder.kmers_present(sequence) {
            return None;
        }
        let alignment = self.core.aligner.locate(sequence)?;
        Some(alignment_to_remove_after(
            alignment,
            sequence,
            &self.core.name,
            &self.core.sequence,
        ))
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn front_params(seq: &str) -> SingleAdapterParams {
        SingleAdapterParams {
            sequence: seq.to_string(),
            max_errors: 0.1,
            min_overlap: 3,
            read_wildcards: false,
            adapter_wildcards: true,
            name: Some("test_front".to_string()),
            indels: true,
        }
    }

    fn back_params(seq: &str) -> SingleAdapterParams {
        SingleAdapterParams {
            sequence: seq.to_string(),
            max_errors: 0.1,
            min_overlap: 3,
            read_wildcards: false,
            adapter_wildcards: true,
            name: Some("test_back".to_string()),
            indels: true,
        }
    }

    #[test]
    fn test_front_adapter_match() {
        let mut adapter = SingleAdapter::Front(FrontAdapter::new(front_params("AGATCGGAAGAGC")));
        assert_eq!(adapter.name(), "test_front");
        assert_eq!(adapter.description(), "regular 5'");
        let m = adapter.match_to("AGATCGGAAGAGCACGTACGT");
        assert!(m.is_some());
        let m = m.unwrap();
        assert!(m.score() > 0);
    }

    #[test]
    fn test_back_adapter_match() {
        let mut adapter = SingleAdapter::Back(BackAdapter::new(back_params("AGATCGGAAGAGC")));
        let m = adapter.match_to("ACGTACGTACGTAGATCGGAAGAGC");
        assert!(m.is_some());
        let m = m.unwrap();
        assert!(m.score() > 0);
    }

    #[test]
    fn test_prefix_adapter_match_no_indels() {
        let params = SingleAdapterParams {
            sequence: "ACGT".to_string(),
            max_errors: 0.0,
            min_overlap: 4,
            read_wildcards: false,
            adapter_wildcards: false,
            name: Some("prefix".to_string()),
            indels: false,
        };
        let mut adapter = SingleAdapter::Prefix(PrefixAdapter::new(params));
        assert_eq!(adapter.description(), "anchored 5'");
        assert!(!adapter.allows_partial_matches());

        let m = adapter.match_to("ACGTTTTT");
        assert!(m.is_some());

        let m = adapter.match_to("TTTTACGT");
        assert!(m.is_none());
    }

    #[test]
    fn test_suffix_adapter_match_no_indels() {
        let params = SingleAdapterParams {
            sequence: "ACGT".to_string(),
            max_errors: 0.0,
            min_overlap: 4,
            read_wildcards: false,
            adapter_wildcards: false,
            name: Some("suffix".to_string()),
            indels: false,
        };
        let mut adapter = SingleAdapter::Suffix(SuffixAdapter::new(params));
        assert_eq!(adapter.description(), "anchored 3'");

        let m = adapter.match_to("TTTTACGT");
        assert!(m.is_some());

        let m = adapter.match_to("ACGTTTTT");
        assert!(m.is_none());
    }

    #[test]
    fn test_anywhere_adapter() {
        let params = SingleAdapterParams {
            sequence: "AGATCGGAAGAGC".to_string(),
            max_errors: 0.1,
            min_overlap: 3,
            read_wildcards: false,
            adapter_wildcards: true,
            name: Some("anywhere".to_string()),
            indels: true,
        };
        let mut adapter = SingleAdapter::Anywhere(AnywhereAdapter::new(params));
        assert_eq!(adapter.description(), "variable 5'/3'");

        // At start -> RemoveBefore
        let m = adapter.match_to("AGATCGGAAGAGCACGTACGT");
        assert!(m.is_some());

        // At end -> RemoveAfter
        let m = adapter.match_to("ACGTACGTACGTAGATCGGAAGAGC");
        assert!(m.is_some());
    }

    #[test]
    fn test_no_match() {
        let mut adapter = SingleAdapter::Front(FrontAdapter::new(front_params("AGATCGGAAGAGC")));
        let m = adapter.match_to("TTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        assert!(m.is_none());
    }

    #[test]
    fn test_adapter_spec() {
        let adapter = SingleAdapter::Front(FrontAdapter::new(front_params("ACGT")));
        assert_eq!(adapter.spec(), "ACGT...");

        let adapter = SingleAdapter::Back(BackAdapter::new(back_params("ACGT")));
        assert_eq!(adapter.spec(), "ACGT");
    }

    #[test]
    fn test_create_statistics() {
        let adapter = SingleAdapter::Front(FrontAdapter::new(front_params("ACGT")));
        let stats = adapter.create_statistics();
        assert_eq!(stats.name(), "test_front");
        let (front, back) = stats.end_statistics();
        assert!(front.is_some());
        assert!(back.is_none());
    }

    #[test]
    fn test_non_internal_front() {
        let params = SingleAdapterParams {
            sequence: "AGATCGGAAGAGC".to_string(),
            max_errors: 0.1,
            min_overlap: 3,
            read_wildcards: false,
            adapter_wildcards: true,
            name: Some("ni_front".to_string()),
            indels: true,
        };
        let adapter =
            SingleAdapter::NonInternalFront(NonInternalFrontAdapter::new(params));
        assert_eq!(adapter.description(), "non-internal 5'");
        assert_eq!(adapter.descriptive_identifier(), "noninternal_five_prime");
    }

    #[test]
    fn test_non_internal_back() {
        let params = SingleAdapterParams {
            sequence: "AGATCGGAAGAGC".to_string(),
            max_errors: 0.1,
            min_overlap: 3,
            read_wildcards: false,
            adapter_wildcards: true,
            name: Some("ni_back".to_string()),
            indels: true,
        };
        let adapter =
            SingleAdapter::NonInternalBack(NonInternalBackAdapter::new(params));
        assert_eq!(adapter.description(), "non-internal 3'");
    }

    #[test]
    fn test_rightmost_front() {
        let params = SingleAdapterParams {
            sequence: "AGATCGGAAGAGC".to_string(),
            max_errors: 0.1,
            min_overlap: 3,
            read_wildcards: false,
            adapter_wildcards: true,
            name: Some("rm_front".to_string()),
            indels: true,
        };
        let mut adapter =
            SingleAdapter::RightmostFront(RightmostFrontAdapter::new(params));
        assert_eq!(adapter.description(), "rightmost 5'");
        let m = adapter.match_to("AGATCGGAAGAGCACGTACGT");
        assert!(m.is_some());
    }

    #[test]
    fn test_rightmost_back() {
        let params = SingleAdapterParams {
            sequence: "AGATCGGAAGAGC".to_string(),
            max_errors: 0.1,
            min_overlap: 3,
            read_wildcards: false,
            adapter_wildcards: true,
            name: Some("rm_back".to_string()),
            indels: true,
        };
        let mut adapter =
            SingleAdapter::RightmostBack(RightmostBackAdapter::new(params));
        assert_eq!(adapter.description(), "rightmost 3'");
        let m = adapter.match_to("ACGTACGTACGTAGATCGGAAGAGC");
        assert!(m.is_some());
    }

    #[test]
    fn test_sequence_normalisation() {
        let params = SingleAdapterParams {
            sequence: "acguacgu".to_string(),
            max_errors: 0.1,
            min_overlap: 3,
            read_wildcards: false,
            adapter_wildcards: true,
            name: Some("norm".to_string()),
            indels: true,
        };
        let adapter = SingleAdapter::Front(FrontAdapter::new(params));
        assert_eq!(adapter.sequence(), "ACGTACGT");
    }
}
