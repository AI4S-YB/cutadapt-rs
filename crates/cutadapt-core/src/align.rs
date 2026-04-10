// Sequence alignment algorithms.
//
// Port of cutadapt/_align.pyx and cutadapt/align.py.

use crate::match_tables::{acgt_table, iupac_table, translate, upper_table};
use bitflags::bitflags;

// Scoring constants
const MATCH_SCORE: i32 = 1;
const MISMATCH_SCORE: i32 = -1;
const INSERTION_SCORE: i32 = -2;
const DELETION_SCORE: i32 = -2;

bitflags! {
    /// Flags for the Aligner that indicate which ends of reference or query
    /// may be skipped at no cost.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct EndSkip: u8 {
        const REFERENCE_START = 1;
        const QUERY_START = 2;
        const REFERENCE_END = 4;
        const QUERY_STOP = 8;
        const SEMIGLOBAL = 15;
    }
}

/// A single cell in the DP matrix column.
#[derive(Debug, Clone, Copy)]
struct Entry {
    cost: i32,
    score: i32,
    origin: i32,
}

/// Best match found during alignment.
#[derive(Debug, Clone, Copy)]
struct Match {
    origin: i32,
    cost: i32,
    score: i32,
    ref_stop: i32,
    query_stop: i32,
}

/// Alignment result tuple: (ref_start, ref_stop, query_start, query_stop, score, errors).
pub type AlignmentResult = (usize, usize, usize, usize, i32, i32);

/// Find a full or partial occurrence of a query string in a reference string
/// allowing errors (mismatches, insertions, deletions).
///
/// Uses a hybrid cost/score semi-global DP alignment algorithm.
pub struct Aligner {
    /// Reference string (as provided by user).
    pub reference: String,
    /// Internal encoded reference bytes.
    encoded_ref: Vec<u8>,
    /// DP matrix column (reused across calls).
    column: Vec<Entry>,
    /// Cumulative N-character counts: n_counts[i] = number of N's in reference[..i].
    n_counts: Vec<i32>,
    /// Effective length (reference length minus N wildcards if wildcard_ref).
    pub effective_length: i32,
    /// Maximum allowed error rate.
    pub max_error_rate: f64,
    /// Alignment flags.
    start_in_reference: bool,
    start_in_query: bool,
    stop_in_reference: bool,
    stop_in_query: bool,
    /// Wildcard options.
    wildcard_ref: bool,
    wildcard_query: bool,
    /// Indel costs.
    insertion_cost: i32,
    deletion_cost: i32,
    /// Minimum overlap requirement.
    min_overlap: i32,
    /// Reference length.
    m: i32,
}

impl Aligner {
    pub fn new(
        reference: &str,
        max_error_rate: f64,
        flags: u8,
        wildcard_ref: bool,
        wildcard_query: bool,
        indel_cost: i32,
        min_overlap: i32,
    ) -> Self {
        assert!(indel_cost >= 1, "indel_cost must be at least 1");
        let m = reference.len() as i32;

        // Compute n_counts
        let mut n_counts = vec![0i32; (m + 1) as usize];
        let mut n_count = 0i32;
        for (i, ch) in reference.chars().enumerate() {
            n_counts[i] = n_count;
            if ch == 'N' || ch == 'n' {
                n_count += 1;
            }
        }
        n_counts[m as usize] = n_count;

        let mut effective_length = m;
        let encoded_ref;
        if wildcard_ref {
            effective_length = m - n_counts[m as usize];
            if effective_length == 0 {
                panic!("Cannot have only N wildcards in the sequence");
            }
            let table = iupac_table();
            encoded_ref = translate(reference.as_bytes(), &table);
        } else if wildcard_query {
            let table = acgt_table();
            encoded_ref = translate(reference.as_bytes(), &table);
        } else {
            encoded_ref = reference.as_bytes().to_ascii_uppercase();
        }

        let column = vec![
            Entry {
                cost: 0,
                score: 0,
                origin: 0,
            };
            (m + 1) as usize
        ];

        Self {
            reference: reference.to_string(),
            encoded_ref,
            column,
            n_counts,
            effective_length,
            max_error_rate,
            start_in_reference: (flags & 1) != 0,
            start_in_query: (flags & 2) != 0,
            stop_in_reference: (flags & 4) != 0,
            stop_in_query: (flags & 8) != 0,
            wildcard_ref,
            wildcard_query,
            insertion_cost: indel_cost,
            deletion_cost: indel_cost,
            min_overlap,
            m,
        }
    }

    /// Find the query within the reference.
    ///
    /// Returns `Some((ref_start, ref_stop, query_start, query_stop, score, errors))`
    /// or `None` if no acceptable alignment was found.
    pub fn locate(&mut self, query: &str) -> Option<AlignmentResult> {
        let m = self.m;
        let n = query.len() as i32;
        let max_error_rate = self.max_error_rate;
        let stop_in_query = self.stop_in_query;

        // Encode query
        let (query_bytes, compare_ascii) = if self.wildcard_query {
            let table = iupac_table();
            (translate(query.as_bytes(), &table), false)
        } else if self.wildcard_ref {
            let table = acgt_table();
            (translate(query.as_bytes(), &table), false)
        } else {
            (query.as_bytes().to_ascii_uppercase(), true)
        };

        let s1 = &self.encoded_ref;
        let s2 = &query_bytes;

        // Maximum number of errors
        let k = (max_error_rate * m as f64) as i32;

        // Determine largest and smallest column we need to compute
        let mut max_n = n;
        let mut min_n = 0i32;
        if !self.start_in_query {
            max_n = n.min(m + k);
        }
        if !self.stop_in_query {
            min_n = 0i32.max(n - m - k);
        }

        // Fill column min_n
        let column = &mut self.column;
        if !self.start_in_reference && !self.start_in_query {
            for i in 0..=(m as usize) {
                let ii = i as i32;
                column[i].score = ii * DELETION_SCORE;
                column[i].cost = ii.max(min_n) * self.deletion_cost;
                column[i].origin = 0;
            }
        } else if self.start_in_reference && !self.start_in_query {
            for i in 0..=(m as usize) {
                let ii = i as i32;
                column[i].score = 0;
                column[i].cost = min_n * self.deletion_cost;
                column[i].origin = (min_n - ii).min(0);
            }
        } else if !self.start_in_reference && self.start_in_query {
            for i in 0..=(m as usize) {
                let ii = i as i32;
                column[i].score = ii * DELETION_SCORE;
                column[i].cost = ii * self.deletion_cost;
                column[i].origin = (min_n - ii).max(0);
            }
        } else {
            // start_in_reference && start_in_query
            for i in 0..=(m as usize) {
                let ii = i as i32;
                column[i].score = 0;
                column[i].cost = ii.min(min_n) * self.deletion_cost;
                column[i].origin = min_n - ii;
            }
        }

        let mut best = Match {
            ref_stop: m,
            query_stop: n,
            cost: m + n + 1,
            origin: 0,
            score: 0,
        };

        // Ukkonen's trick: index of the last cell that is at most k
        let mut last = if self.start_in_reference {
            m
        } else {
            m.min(k + 1)
        };

        let origin_increment: i32 = if self.start_in_query { 1 } else { 0 };
        let insertion_cost_increment: i32 = if self.start_in_query {
            0
        } else {
            self.insertion_cost
        };
        let insertion_score_increment: i32 = if self.start_in_query { 0 } else { INSERTION_SCORE };

        let insertion_cost = self.insertion_cost;
        let deletion_cost = self.deletion_cost;
        let match_score = MATCH_SCORE;
        let mismatch_score = MISMATCH_SCORE;
        let insertion_score = INSERTION_SCORE;
        let deletion_score = DELETION_SCORE;

        let mut last_filled_i: i32 = 0;

        // Iterate over columns
        for j in (min_n + 1)..=(max_n) {
            let ju = j as usize;

            // Remember first entry before overwriting
            let mut diag_entry = column[0];

            // Fill in first entry in this column
            column[0].origin += origin_increment;
            column[0].cost += insertion_cost_increment;
            column[0].score += insertion_score_increment;

            for i in 1..=(last as usize) {
                let characters_equal = if compare_ascii {
                    s1[i - 1] == s2[ju - 1]
                } else {
                    (s1[i - 1] & s2[ju - 1]) != 0
                };

                let (cost, origin, score) = if characters_equal {
                    (
                        diag_entry.cost,
                        diag_entry.origin,
                        diag_entry.score + match_score,
                    )
                } else {
                    let current_entry = column[i];
                    let previous_entry = column[i - 1];
                    let cost_diag = diag_entry.cost + 1;
                    let cost_insertion = current_entry.cost + insertion_cost;
                    let cost_deletion = previous_entry.cost + deletion_cost;

                    if cost_diag <= cost_deletion && cost_diag <= cost_insertion {
                        // MISMATCH
                        (cost_diag, diag_entry.origin, diag_entry.score + mismatch_score)
                    } else if cost_deletion <= cost_insertion {
                        // DELETION
                        (
                            cost_deletion,
                            previous_entry.origin,
                            previous_entry.score + deletion_score,
                        )
                    } else {
                        // INSERTION
                        (
                            cost_insertion,
                            current_entry.origin,
                            current_entry.score + insertion_score,
                        )
                    }
                };

                // Remember the current cell for next iteration
                diag_entry = column[i];

                column[i].cost = cost;
                column[i].origin = origin;
                column[i].score = score;
            }
            last_filled_i = last;

            while last >= 0 && column[last as usize].cost > k {
                last -= 1;
            }
            if last < m {
                last += 1;
            } else if stop_in_query {
                // Found a match. Find best match in last row.
                let cost = column[m as usize].cost;
                let score = column[m as usize].score;
                let origin = column[m as usize].origin;
                let length = m + origin.min(0);
                let cur_effective_length = if self.wildcard_ref {
                    if length < m {
                        length
                            - (self.n_counts[m as usize]
                                - self.n_counts[(m - length) as usize])
                    } else {
                        self.effective_length
                    }
                } else {
                    length
                };

                let is_acceptable = length >= self.min_overlap
                    && cost as f64 <= cur_effective_length as f64 * max_error_rate;

                let best_length = m + best.origin.min(0);

                if is_acceptable
                    && ((best.cost == m + n + 1)
                        || (origin <= best.origin + m / 2 && score > best.score)
                        || (length > best_length && score > best.score))
                {
                    best.score = score;
                    best.cost = cost;
                    best.origin = origin;
                    best.ref_stop = m;
                    best.query_stop = j;
                    if cost == 0 && origin >= 0 {
                        // Exact match, stop early
                        break;
                    }
                }
            }
        }

        if max_n == n {
            let first_i: i32 = if self.stop_in_reference { 0 } else { m };
            // Search in last column (in reverse)
            for i in (first_i..=(last_filled_i)).rev() {
                let iu = i as usize;
                let length = i + column[iu].origin.min(0);
                let cost = column[iu].cost;
                let score = column[iu].score;
                let cur_effective_length = if self.wildcard_ref {
                    if length < m {
                        let ref_start = (-column[iu].origin).max(0);
                        length - (self.n_counts[iu] - self.n_counts[ref_start as usize])
                    } else {
                        self.effective_length
                    }
                } else {
                    length
                };

                let is_acceptable = length >= self.min_overlap
                    && cost as f64 <= cur_effective_length as f64 * max_error_rate;

                let best_length = best.ref_stop + best.origin.min(0);

                if is_acceptable
                    && ((best.cost == m + n + 1)
                        || (column[iu].origin <= best.origin + m / 2 && score > best.score)
                        || (length > best_length && score > best.score))
                {
                    best.score = score;
                    best.cost = cost;
                    best.origin = column[iu].origin;
                    best.ref_stop = i;
                    best.query_stop = n;
                }
            }
        }

        if best.cost == m + n + 1 {
            return None;
        }

        let (ref_start, query_start) = if best.origin >= 0 {
            (0usize, best.origin as usize)
        } else {
            ((-best.origin) as usize, 0usize)
        };

        Some((
            ref_start,
            best.ref_stop as usize,
            query_start,
            best.query_stop as usize,
            best.score,
            best.cost,
        ))
    }
}

/// Fast prefix comparison (no indels, Hamming distance only).
/// Used for anchored 5' adapters in no-indels mode.
pub struct PrefixComparer {
    encoded_ref: Vec<u8>,
    wildcard_ref: bool,
    wildcard_query: bool,
    m: usize,
    max_k: i32,
    pub effective_length: i32,
    min_overlap: usize,
}

impl PrefixComparer {
    pub fn new(
        reference: &str,
        max_error_rate: f64,
        wildcard_ref: bool,
        wildcard_query: bool,
        min_overlap: usize,
    ) -> Self {
        let m = reference.len();
        let mut effective_length = m as i32;
        if wildcard_ref {
            effective_length -= reference.matches('N').count() as i32
                + reference.matches('n').count() as i32;
            if effective_length == 0 {
                panic!("Cannot have only N wildcards in the sequence");
            }
        }
        let max_k = (max_error_rate * effective_length as f64) as i32;

        let encoded_ref = if wildcard_ref {
            let table = iupac_table();
            translate(reference.as_bytes(), &table)
        } else if wildcard_query {
            let table = acgt_table();
            translate(reference.as_bytes(), &table)
        } else {
            let table = upper_table();
            translate(reference.as_bytes(), &table)
        };

        Self {
            encoded_ref,
            wildcard_ref,
            wildcard_query,
            m,
            max_k,
            effective_length,
            min_overlap,
        }
    }

    pub fn locate(&self, query: &str) -> Option<AlignmentResult> {
        let n = query.len();
        let length = self.m.min(n);

        let (query_bytes, compare_ascii) = if self.wildcard_query {
            let table = iupac_table();
            (translate(query.as_bytes(), &table), false)
        } else if self.wildcard_ref {
            let table = acgt_table();
            (translate(query.as_bytes(), &table), false)
        } else {
            let table = upper_table();
            (translate(query.as_bytes(), &table), true)
        };

        let mut errors = 0i32;
        if compare_ascii {
            for i in 0..length {
                if self.encoded_ref[i] != query_bytes[i] {
                    errors += 1;
                }
            }
        } else {
            for i in 0..length {
                if (self.encoded_ref[i] & query_bytes[i]) == 0 {
                    errors += 1;
                }
            }
        }

        if errors > self.max_k || length < self.min_overlap {
            return None;
        }
        let score = (length as i32 - errors) * MATCH_SCORE + errors * MISMATCH_SCORE;
        Some((0, length, 0, length, score, errors))
    }
}

/// Fast suffix comparison (no indels, Hamming distance only).
/// Used for anchored 3' adapters in no-indels mode.
pub struct SuffixComparer {
    inner: PrefixComparer,
    m: usize,
}

impl SuffixComparer {
    pub fn new(
        reference: &str,
        max_error_rate: f64,
        wildcard_ref: bool,
        wildcard_query: bool,
        min_overlap: usize,
    ) -> Self {
        let reversed: String = reference.chars().rev().collect();
        let m = reference.len();
        Self {
            inner: PrefixComparer::new(&reversed, max_error_rate, wildcard_ref, wildcard_query, min_overlap),
            m,
        }
    }

    pub fn locate(&self, query: &str) -> Option<AlignmentResult> {
        let n = query.len();
        let reversed_query: String = query.chars().rev().collect();
        let result = self.inner.locate(&reversed_query)?;
        let (_, length, _, _, score, errors) = result;
        Some((self.m - length, self.m, n - length, n, score, errors))
    }

    pub fn effective_length(&self) -> i32 {
        self.inner.effective_length
    }
}

/// Yield all strings t for which the Hamming distance between s and t is exactly k,
/// assuming the alphabet is A, C, G, T.
pub fn hamming_sphere(s: &str, k: usize) -> Vec<String> {
    let mut results = Vec::new();
    hamming_sphere_recursive(s.as_bytes(), k, &mut results);
    results
}

fn hamming_sphere_recursive(s: &[u8], k: usize, results: &mut Vec<String>) {
    if k == 0 {
        results.push(String::from_utf8_lossy(s).to_string());
        return;
    }
    let n = s.len();
    let alphabet = b"ACGT";

    if k == 1 {
        for i in 0..n {
            for &ch in alphabet {
                if s[i] == ch {
                    continue;
                }
                let mut result = s.to_vec();
                result[i] = ch;
                results.push(String::from_utf8_lossy(&result).to_string());
            }
        }
        return;
    }

    if k == 2 {
        for i in 0..n {
            for &ch1 in alphabet {
                if s[i] == ch1 {
                    continue;
                }
                for j in (i + 1)..n {
                    for &ch2 in alphabet {
                        if s[j] == ch2 {
                            continue;
                        }
                        let mut result = s.to_vec();
                        result[i] = ch1;
                        result[j] = ch2;
                        results.push(String::from_utf8_lossy(&result).to_string());
                    }
                }
            }
        }
        return;
    }

    // Recursive for k > 2
    for i in 0..=(n - k) {
        let prefix = &s[..i];
        let c = s[i];
        let suffix = &s[i + 1..];
        for &pch in alphabet {
            if pch == c {
                continue;
            }
            let mut sub_results = Vec::new();
            hamming_sphere_recursive(suffix, k - 1, &mut sub_results);
            for t in sub_results {
                let mut full = Vec::with_capacity(n);
                full.extend_from_slice(prefix);
                full.push(pch);
                full.extend_from_slice(t.as_bytes());
                results.push(String::from_utf8_lossy(&full).to_string());
            }
        }
    }
}

/// Yield all strings t for which the Hamming distance between s and t is at most k.
/// Returns tuples (t, distance, matches) where matches = len(t) - distance.
pub fn hamming_environment(s: &str, k: usize) -> Vec<(String, usize, usize)> {
    let n = s.len();
    let mut results = Vec::new();
    for e in 0..=k {
        for t in hamming_sphere(s, e) {
            results.push((t, e, n - e));
        }
    }
    results
}

/// Return the edit distance between strings s and t.
pub fn edit_distance(s: &str, t: &str) -> usize {
    let s = s.as_bytes();
    let t = t.as_bytes();
    let m = s.len();
    let n = t.len();
    let mut costs: Vec<usize> = (0..=m).collect();

    for j in 1..=n {
        let mut prev = costs[0];
        costs[0] += 1;
        for i in 1..=m {
            let match_cost = if s[i - 1] == t[j - 1] { 0 } else { 1 };
            let c = (prev + match_cost).min(costs[i] + 1).min(costs[i - 1] + 1);
            prev = costs[i];
            costs[i] = c;
        }
    }
    costs[m]
}

/// Find all strings s for which the edit distance between s and t is at most k,
/// assuming the alphabet is A, C, G, T.
///
/// Returns Vec of (s, edit_distance, num_matches).
pub fn edit_environment(t: &str, k: usize) -> Vec<(String, usize, usize)> {
    let t_bytes: Vec<u8> = t
        .bytes()
        .map(|b| match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => b,
        })
        .collect();

    let n = t_bytes.len();
    let rows = n + k + 1;
    let cols = n + 1;

    // costs[i * cols + j]
    let mut costs = vec![(k + 1) as i32; rows * cols];
    for i in 0..rows {
        costs[i * cols] = i as i32;
    }
    for j in 0..cols {
        costs[j] = j as i32;
    }

    let mut matches = vec![0i32; rows * cols];

    let trans = [b'A', b'C', b'G', b'T'];
    let mut results = Vec::new();

    let mut s = vec![0u8; rows];
    let mut i: usize = 0;

    loop {
        // Fill in row i
        if i > 0 {
            let ch = s[i - 1];
            let j_start = if i > k { i - k } else { 1 };
            let j_end = (n + 1).min(i + k + 1);
            for j in j_start..j_end {
                let match_cost = if t_bytes[j - 1] == ch { 0i32 } else { 1i32 };
                let diag = costs[(i - 1) * cols + j - 1] + match_cost;
                let left = costs[i * cols + j - 1] + 1;
                let up = costs[(i - 1) * cols + j] + 1;

                let (c, m) = if diag <= left && diag <= up {
                    (diag, matches[(i - 1) * cols + j - 1] + (1 - match_cost))
                } else if left <= up {
                    (left, matches[i * cols + j - 1])
                } else {
                    (up, matches[(i - 1) * cols + j])
                };
                costs[i * cols + j] = c;
                matches[i * cols + j] = m;
            }
        }

        let min_cost = if i > 0 {
            let j_start = if i > k { i - k } else { 1 };
            let j_end = (n + 1).min(i + k + 1);
            (j_start..j_end)
                .map(|j| costs[i * cols + j])
                .min()
                .unwrap_or(k as i32 + 1)
        } else {
            0
        };

        if costs[i * cols + n] <= k as i32 {
            let result_str: String = s[..i].iter().map(|&b| trans[b as usize] as char).collect();
            results.push((
                result_str,
                costs[i * cols + n] as usize,
                matches[i * cols + n] as usize,
            ));
        }

        // Next string
        if min_cost <= k as i32 && i < n + k {
            s[i] = 0;
            i += 1;
        } else {
            loop {
                if i == 0 {
                    return results;
                }
                i -= 1;
                let ch = s[i];
                if ch < 3 {
                    s[i] = ch + 1;
                    i += 1;
                    break;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn locate(
        reference: &str,
        query: &str,
        max_error_rate: f64,
        flags: u8,
        wildcard_ref: bool,
        wildcard_query: bool,
        min_overlap: i32,
    ) -> Option<AlignmentResult> {
        let mut aligner = Aligner::new(
            reference,
            max_error_rate,
            flags,
            wildcard_ref,
            wildcard_query,
            1,
            min_overlap,
        );
        aligner.locate(query)
    }

    fn locate_default(
        reference: &str,
        query: &str,
        max_error_rate: f64,
        flags: u8,
    ) -> Option<AlignmentResult> {
        locate(reference, query, max_error_rate, flags, false, false, 1)
    }

    // Where flags (from adapters.py Where enum)
    const WHERE_BACK: u8 = EndSkip::QUERY_START.bits() | EndSkip::QUERY_STOP.bits() | EndSkip::REFERENCE_END.bits(); // 14
    const WHERE_FRONT: u8 = EndSkip::QUERY_START.bits() | EndSkip::QUERY_STOP.bits() | EndSkip::REFERENCE_START.bits(); // 11
    const WHERE_PREFIX: u8 = EndSkip::QUERY_STOP.bits(); // 8

    #[test]
    fn test_basic() {
        let reference = "CTCCAGCTTAGACATATC";
        let mut aligner = Aligner::new(reference, 0.1, WHERE_BACK, false, false, 1, 1);
        aligner.locate("CC");
    }

    #[test]
    fn test_100_percent_error_rate() {
        let reference = "GCTTAGACATATC";
        let mut aligner = Aligner::new(reference, 1.0, WHERE_BACK, false, false, 1, 1);
        aligner.locate("CAA");
    }

    #[test]
    #[should_panic(expected = "Cannot have only N wildcards")]
    fn test_not_only_n_wildcards() {
        Aligner::new("NNNNN", 0.1, 15, true, false, 1, 1);
    }

    #[test]
    fn test_find_empty_in_empty() {
        let mut aligner = Aligner::new("", 0.0, 0, false, false, 1, 0);
        let result = aligner.locate("");
        assert_eq!(result, Some((0, 0, 0, 0, 0, 0)));
    }

    #[test]
    fn test_indels_penalized() {
        let mut aligner = Aligner::new("CCAGTCCTCT", 0.3, WHERE_PREFIX, false, false, 1, 1);
        let result = aligner.locate("CCAGTCCTTTCCTGAGAGT");
        assert_eq!(result, Some((0, 10, 0, 10, 8, 1)));

        let mut aligner = Aligner::new("TCGATC", (1.5 / 6.0), WHERE_PREFIX, false, false, 1, 1);
        let result = aligner.locate("TCGATGC");
        assert_eq!(result, Some((0, 6, 0, 6, 4, 1)));
    }

    #[test]
    fn test_align_illumina() {
        let mut aligner =
            Aligner::new("GCCGAACTTCTTAGACTGCCTTAAGGACGT", 0.1, WHERE_BACK, false, false, 1, 1);
        let result = aligner
            .locate("CAAATCACCAGAAGGCGCCTAACTTCTTAGACTGCC")
            .unwrap();
        assert_eq!(result.0, 0); // ref_start
        assert_eq!(result.1, 20); // ref_end
        assert_eq!(result.2, 16); // query_start
        // query_end: the Python test says 36, but query length is 35
        // Let's check: "CAAATCACCAGAAGGCGCCTAACTTCTTAGACTGCC" has 35 chars
        assert_eq!(result.3, 36); // query_end
        assert_eq!(result.4, 18); // score
        assert_eq!(result.5, 1); // errors
    }

    #[test]
    fn test_poly_t() {
        let mut aligner = Aligner::new("TTTT", 0.25, WHERE_BACK, false, false, 1, 1);
        let result = aligner.locate("CCTTTT").unwrap();
        assert_eq!(result, (0, 4, 2, 6, 4, 0));
    }

    #[test]
    fn test_poly_t_partial_match() {
        let mut aligner = Aligner::new("TTTTTT", 0.25, WHERE_BACK, false, false, 1, 1);
        let result = aligner.locate("CCTTTT").unwrap();
        assert_eq!(result, (0, 4, 2, 6, 4, 0));
    }

    #[test]
    fn test_poly_t_2() {
        let mut aligner = Aligner::new("TTT", 1.0 / 3.0, WHERE_BACK, false, false, 1, 1);
        let result = aligner.locate("CCTTTT").unwrap();
        assert_eq!(result.0, 0); // ref_start
        assert_eq!(result.1, 3); // ref_stop
        assert_eq!(result.2, 2); // query_start
        assert_eq!(result.3, 5); // query_stop
    }

    #[test]
    fn test_poly_a() {
        let s = "AAAAAAAAAAAAAAAAA";
        let t = "ACAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let result = locate_default(s, t, 0.0, WHERE_BACK).unwrap();
        assert_eq!(
            result,
            (0, s.len(), 4, 4 + s.len(), s.len() as i32, 0)
        );
    }

    #[test]
    fn test_no_match() {
        let result = locate_default("CTGATCTGGCCG", "AAAAGGG", 0.1, WHERE_BACK);
        assert!(result.is_none());
    }

    #[test]
    fn test_wildcards_in_adapter() {
        let wildcard_sequences = vec![
            "CCCATTGATC",
            "CCCRTTRATC",
            "YCCATYGATC",
            "CSSATTSATC",
            "CCCWWWGATC",
            "CCCATKKATC",
            "CCMATTGMTC",
            "BCCATTBABC",
            "BCCATTBABC",
            "CCCDTTDADC",
            "CHCATHGATC",
            "CVCVTTVATC",
            "CCNATNGATC",
            "CCCNTTNATC",
        ];
        let original = wildcard_sequences[0];
        let r = format!("CATCTGTCC{}GCCAGGGTTGATTCGGCTGATCTGGCCG", original);
        for a in &wildcard_sequences {
            let result = locate(a, &r, 0.0, WHERE_BACK, true, false, 1).unwrap();
            assert_eq!(result, (0, 10, 9, 19, 10, 0));
        }
        // X should not match
        let result = locate("CCCXTTXATC", &r, 0.0, WHERE_BACK, true, false, 1);
        assert!(result.is_none());
    }

    #[test]
    fn test_wildcards_in_read() {
        let wildcard_sequences = vec![
            "CCCATTGATC",
            "CCCRTTRATC",
            "YCCATYGATC",
            "CSSATTSATC",
            "CCCWWWGATC",
            "CCCATKKATC",
            "CCMATTGMTC",
            "BCCATTBABC",
            "BCCATTBABC",
            "CCCDTTDADC",
            "CHCATHGATC",
            "CVCVTTVATC",
            "CCNATNGATC",
            "CCCNTTNATC",
            "CCCXTTXATC",
        ];
        let a = "CCCATTGATC";
        for s in &wildcard_sequences {
            let r = format!("CATCTGTCC{}GCCAGGGTTGATTCGGCTGATCTGGCCG", s);
            let result = locate(a, &r, 0.0, WHERE_BACK, false, true, 1);
            if s.contains('X') {
                assert!(result.is_none());
            } else {
                assert_eq!(result.unwrap(), (0, 10, 9, 19, 10, 0));
            }
        }
    }

    #[test]
    fn test_wildcards_in_both() {
        let wildcard_sequences = vec![
            "CCCATTGATC",
            "CCCRTTRATC",
            "YCCATYGATC",
            "CSSATTSATC",
            "CCCWWWGATC",
            "CCCATKKATC",
            "CCMATTGMTC",
            "BCCATTBABC",
            "BCCATTBABC",
            "CCCDTTDADC",
            "CHCATHGATC",
            "CVCVTTVATC",
            "CCNATNGATC",
            "CCCNTTNATC",
        ];
        for a in &wildcard_sequences {
            for s in &wildcard_sequences {
                let r = format!("CATCTGTCC{}GCCAGGGTTGATTCGGCTGATCTGGCCG", s);
                let result = locate(a, &r, 0.0, WHERE_BACK, true, true, 1).unwrap();
                assert_eq!(result, (0, 10, 9, 19, 10, 0));
            }
        }
    }

    #[test]
    fn test_n_wildcard_in_ref_matches_n_in_query_prefix() {
        let comparer = PrefixComparer::new("NNACGT", 0.0, true, false, 1);
        let result = comparer.locate("NTACGTAA");
        assert_eq!(result, Some((0, 6, 0, 6, 6, 0)));

        let result2 = comparer.locate("YTACGTAA");
        assert_eq!(result2, Some((0, 6, 0, 6, 6, 0)));
    }

    #[test]
    fn test_n_wildcard_in_ref_matches_n_in_query_back() {
        let mut aligner = Aligner::new("NNACGT", 0.0, WHERE_BACK, true, false, 1, 1);
        let result = aligner.locate("AAANTACGTAAA").unwrap();
        assert_eq!(result, (0, 6, 3, 9, 6, 0));
    }

    #[test]
    fn test_compare_prefixes() {
        let c = PrefixComparer::new("AAXAA", 0.9, false, false, 1);
        assert_eq!(c.locate("AAAAATTTTTTTTT"), Some((0, 5, 0, 5, 3, 1)));

        let c = PrefixComparer::new("AANAA", 0.9, true, false, 1);
        assert_eq!(c.locate("AACAATTTTTTTTT"), Some((0, 5, 0, 5, 5, 0)));

        let c = PrefixComparer::new("XAAAAA", 0.9, false, false, 1);
        assert_eq!(c.locate("AAAAATTTTTTTTT"), Some((0, 6, 0, 6, 2, 2)));
    }

    #[test]
    fn test_compare_suffixes() {
        let c = SuffixComparer::new("AAXAA", 0.9, false, false, 1);
        assert_eq!(c.locate("TTTTTTTAAAAA"), Some((0, 5, 7, 12, 3, 1)));

        let c = SuffixComparer::new("AANAA", 0.9, true, false, 1);
        assert_eq!(c.locate("TTTTTTTAACAA"), Some((0, 5, 7, 12, 5, 0)));

        let c = SuffixComparer::new("AAAAAX", 0.9, false, false, 1);
        assert_eq!(c.locate("TTTTTTTAAAAA"), Some((0, 6, 6, 12, 2, 2)));
    }

    #[test]
    fn test_prefix_comparer_none() {
        let c = PrefixComparer::new("AXCGT", 0.4, false, false, 1);
        assert!(c.locate("TTG").is_none());
        assert!(c.locate("AGT").is_some());
        assert!(c.locate("agt").is_some());
        assert!(c.locate("CGT").is_none());
    }

    #[test]
    fn test_suffix_comparer_none() {
        let c = SuffixComparer::new("AXCGT", 0.4, false, false, 1);
        assert!(c.locate("TTG").is_none());
        assert!(c.locate("AGT").is_some());
        assert!(c.locate("agt").is_some());
        assert!(c.locate("CGT").is_some());
        assert!(c.locate("TTG").is_none());
    }

    #[test]
    fn test_n_wildcards_not_counted_prefix() {
        let c = PrefixComparer::new("CNNNNNNNNGTT", 0.25, true, false, 1);
        assert!(c.locate("CAAAAAAAAGTT").is_some());
        assert!(c.locate("CAAAAAAAAGTA").is_some());
        assert!(c.locate("CAAAAAAAAGAA").is_none());
    }

    #[test]
    fn test_n_wildcards_not_counted_suffix() {
        let c = SuffixComparer::new("CNNNNNNNNGTT", 0.25, true, false, 1);
        assert!(c.locate("CAAAAAAAAGTT").is_some());
        assert!(c.locate("CAAAAAAAAGTA").is_some());
        assert!(c.locate("CAAAAAAAAGAA").is_none());
    }

    #[test]
    fn test_n_wildcards_not_counted_aligner_back() {
        let ref_seq = "AGGNNNNNNNNNNNNNNTTC";
        let mut aligner = Aligner::new(ref_seq, 0.1, WHERE_BACK, true, false, 1, 3);
        assert_eq!(aligner.effective_length, 6);
        assert!(aligner.locate("TTC").is_none());
        let result = aligner.locate("AGG").unwrap();
        assert_eq!((result.0, result.1, result.2, result.3), (0, 3, 0, 3));
        let result = aligner.locate("AGGCCCCCCC").unwrap();
        assert_eq!((result.0, result.1, result.2, result.3), (0, 10, 0, 10));
        assert!(aligner.locate("ATGCCCCCCC").is_none());
        assert!(aligner.locate("AGGCCCCCCCCCCCCCCATC").is_none());
        let full_ref = ref_seq.replace('N', "G");
        let query = format!("CCC{}AAA", full_ref);
        let result = aligner.locate(&query).unwrap();
        assert_eq!(result, (0, 20, 3, 23, 20, 0));
    }

    #[test]
    fn test_n_wildcards_not_counted_aligner_front() {
        let ref_seq = "AGGNNNNNNNNNNNNNNTTC";
        let mut aligner = Aligner::new(ref_seq, 0.1, WHERE_FRONT, true, false, 1, 3);
        assert_eq!(aligner.effective_length, 6);
        let result = aligner.locate("TTC").unwrap();
        assert_eq!((result.0, result.1, result.2, result.3), (17, 20, 0, 3));
        assert!(aligner.locate("TGC").is_none());
        let result = aligner.locate("CCCCCCCTTC").unwrap();
        assert_eq!((result.0, result.1, result.2, result.3), (10, 20, 0, 10));
        assert!(aligner.locate("CCCCCCCGTC").is_none());
        let full_ref = ref_seq.replace('N', "G");
        let query = format!("CCC{}AAA", full_ref);
        let result = aligner.locate(&query).unwrap();
        assert_eq!(result, (0, 20, 3, 23, 20, 0));
    }

    #[test]
    fn test_hamming_sphere_explicit() {
        assert_eq!(hamming_sphere("", 0), vec![""]);
        assert_eq!(hamming_sphere("A", 0), vec!["A"]);
        assert_eq!(hamming_sphere("A", 1), vec!["C", "G", "T"]);
        assert_eq!(hamming_sphere("GTC", 0), vec!["GTC"]);
        assert_eq!(
            hamming_sphere("GTC", 1),
            vec!["ATC", "CTC", "TTC", "GAC", "GCC", "GGC", "GTA", "GTG", "GTT"]
        );
    }

    fn binomial(n: usize, k: usize) -> usize {
        if k > n {
            return 0;
        }
        if k == 0 || k == n {
            return 1;
        }
        let k = k.min(n - k);
        let mut result = 1;
        for i in 0..k {
            result = result * (n - i) / (i + 1);
        }
        result
    }

    fn hamming_distance(s: &str, t: &str) -> usize {
        s.chars().zip(t.chars()).filter(|(a, b)| a != b).count()
    }

    #[test]
    fn test_hamming_sphere_properties() {
        let cases = vec![
            ("", 0),
            ("A", 0),
            ("AAA", 1),
            ("ACC", 2),
            ("TCATTA", 3),
            ("AAAAAAAA", 1),
        ];
        for (s, k) in cases {
            let result = hamming_sphere(s, k);
            let result_set: std::collections::HashSet<_> = result.iter().cloned().collect();
            assert_eq!(result.len(), result_set.len(), "Duplicates for ({}, {})", s, k);
            let expected_count = 3usize.pow(k as u32) * binomial(s.len(), k);
            assert_eq!(result.len(), expected_count, "Wrong count for ({}, {})", s, k);
            for t in &result {
                assert_eq!(hamming_distance(s, t), k, "Wrong distance for ({}, {})", s, k);
            }
        }
    }

    #[test]
    fn test_edit_distance() {
        assert_eq!(edit_distance("", ""), 0);
        assert_eq!(edit_distance("A", "A"), 0);
        assert_eq!(edit_distance("A", ""), 1);
        assert_eq!(edit_distance("", "A"), 1);
        assert_eq!(edit_distance("A", "C"), 1);
        assert_eq!(edit_distance("ACGT", "ACGT"), 0);
        assert_eq!(edit_distance("ACGT", "ACGG"), 1);
        assert_eq!(edit_distance("ACGT", "ACG"), 1);
        assert_eq!(edit_distance("ACG", "ACGT"), 1);
        assert_eq!(edit_distance("KITTEN", "SITTING"), 3);
    }

    #[test]
    fn test_edit_environment_basic() {
        let results = edit_environment("A", 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0], ("A".to_string(), 0, 1));

        let results = edit_environment("A", 1);
        let strings: Vec<&str> = results.iter().map(|(s, _, _)| s.as_str()).collect();
        assert!(strings.contains(&"A"));
        assert!(strings.contains(&"C"));
        assert!(strings.contains(&"G"));
        assert!(strings.contains(&"T"));
        assert!(strings.contains(&"")); // deletion
        // Insertions: AA, CA, GA, TA, AC, CC, GC, TC, AG, CG, GG, TG, AT, CT, GT, TT
    }

    #[test]
    fn test_edit_environment_no_duplicates() {
        let cases = vec![
            (0, ""),
            (0, "A"),
            (1, "AAA"),
            (1, "TCATTAGA"),
            (2, "ACC"),
            (3, "TCATTA"),
        ];
        for (k, s) in cases {
            let results = edit_environment(s, k);
            let strings: Vec<&str> = results.iter().map(|(s, _, _)| s.as_str()).collect();
            let set: std::collections::HashSet<&str> = strings.iter().copied().collect();
            assert_eq!(strings.len(), set.len(), "Duplicates for ({}, {})", s, k);
        }
    }
}
