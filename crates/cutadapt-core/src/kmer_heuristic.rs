// K-mer set generation for adapter pre-screening.
//
// Port of cutadapt/kmer_heuristic.py.

use std::collections::{HashMap, HashSet};

/// Partition a sequence into almost equal sized chunks.
/// Returns the set of unique chunks.
pub fn kmer_chunks(sequence: &str, chunks: usize) -> HashSet<String> {
    let n = sequence.len();
    let chunk_size = n / chunks;
    let remainder = n % chunks;

    let mut chunk_sizes: Vec<usize> = Vec::with_capacity(chunks);
    for _ in 0..remainder {
        chunk_sizes.push(chunk_size + 1);
    }
    for _ in 0..(chunks - remainder) {
        chunk_sizes.push(chunk_size);
    }

    let mut offset = 0;
    let mut chunk_set = HashSet::new();
    for size in chunk_sizes {
        chunk_set.insert(sequence[offset..offset + size].to_string());
        offset += size;
    }
    chunk_set
}

/// A SearchSet is (start, stop, set of kmers).
type SearchSet = (isize, Option<isize>, HashSet<String>);

/// Minimize the kmer search list by merging duplicate kmers with different positions.
fn minimize_kmer_search_list(
    kmer_search_list: &[(String, isize, Option<isize>)],
) -> Vec<(String, isize, Option<isize>)> {
    let mut kmer_and_offsets: HashMap<String, Vec<(isize, Option<isize>)>> = HashMap::new();
    for (kmer, start, stop) in kmer_search_list {
        kmer_and_offsets
            .entry(kmer.clone())
            .or_default()
            .push((*start, *stop));
    }

    let mut result = Vec::new();
    for (kmer, positions) in &kmer_and_offsets {
        if positions.len() == 1 {
            let (start, stop) = positions[0];
            result.push((kmer.clone(), start, stop));
            continue;
        }
        if positions.contains(&(0, None)) {
            result.push((kmer.clone(), 0, None));
            continue;
        }
        let front_searches: Vec<_> = positions.iter().filter(|(s, _)| *s == 0).collect();
        let back_searches: Vec<_> = positions.iter().filter(|(_, stop)| stop.is_none()).collect();

        if !front_searches.is_empty() {
            let max_stop = front_searches
                .iter()
                .filter_map(|(_, stop)| *stop)
                .max()
                .unwrap();
            result.push((kmer.clone(), 0, Some(max_stop)));
        }
        if !back_searches.is_empty() {
            let min_start = back_searches.iter().map(|(start, _)| *start).min().unwrap();
            result.push((kmer.clone(), min_start, None));
        }
    }
    result
}

/// Remove kmers that are searched in multiple search sets and merge them.
fn remove_redundant_kmers(
    search_sets: &[SearchSet],
) -> Vec<(isize, Option<isize>, Vec<String>)> {
    let mut kmer_search_list = Vec::new();
    for (start, stop, kmer_set) in search_sets {
        for kmer in kmer_set {
            kmer_search_list.push((kmer.clone(), *start, *stop));
        }
    }
    let minimized = minimize_kmer_search_list(&kmer_search_list);
    let mut result_dict: HashMap<(isize, Option<isize>), Vec<String>> = HashMap::new();
    for (kmer, start, stop) in minimized {
        result_dict.entry((start, stop)).or_default().push(kmer);
    }
    result_dict
        .into_iter()
        .map(|((start, stop), kmers)| (start, stop, kmers))
        .collect()
}

/// Create back-overlap search sets for an adapter.
fn create_back_overlap_searchsets(
    adapter: &str,
    min_overlap: usize,
    error_rate: f64,
) -> Vec<SearchSet> {
    let adapter_length = adapter.len();
    let mut error_lengths: Vec<(usize, usize)> = Vec::new();
    let mut max_error = 0usize;
    let mut search_sets: Vec<SearchSet> = Vec::new();

    for i in 0..=adapter_length {
        if (i as f64 * error_rate) as usize > max_error {
            error_lengths.push((max_error, i - 1));
            max_error += 1;
        }
    }
    error_lengths.push((max_error, adapter_length));

    let mut minimum_length = min_overlap;
    for &(max_errors, length) in &error_lengths {
        if minimum_length > length {
            continue;
        }
        if max_errors == 0 {
            let min_overlap_kmer_length = 5.min(adapter_length);
            if minimum_length < min_overlap_kmer_length {
                for i in minimum_length..min_overlap_kmer_length {
                    let mut kmer_set = HashSet::new();
                    kmer_set.insert(adapter[..i].to_string());
                    search_sets.push((-(i as isize), None, kmer_set));
                }
                minimum_length = min_overlap_kmer_length;
            }
        }
        let kmer_sets = kmer_chunks(&adapter[..minimum_length], max_errors + 1);
        search_sets.push((-(length as isize), None, kmer_sets));
        minimum_length = length + 1;
    }
    search_sets
}

/// Create a set of position and kmer combinations for adapter pre-screening.
///
/// Returns `(start, stop, kmers)` tuples where at least one kmer must be present
/// at the specified position range for an alignment to succeed.
pub fn create_positions_and_kmers(
    adapter: &str,
    min_overlap: usize,
    error_rate: f64,
    back_adapter: bool,
    front_adapter: bool,
    internal: bool,
) -> Vec<(isize, Option<isize>, Vec<String>)> {
    let max_errors = (adapter.len() as f64 * error_rate) as usize;
    let mut search_sets: Vec<SearchSet> = Vec::new();

    if back_adapter {
        search_sets.extend(create_back_overlap_searchsets(
            adapter,
            min_overlap,
            error_rate,
        ));
    }
    if front_adapter {
        let reversed: String = adapter.chars().rev().collect();
        let reversed_back = create_back_overlap_searchsets(&reversed, min_overlap, error_rate);
        for (start, _stop, kmer_set) in reversed_back {
            let new_kmer_set: HashSet<String> = kmer_set
                .iter()
                .map(|k| k.chars().rev().collect())
                .collect();
            search_sets.push((0, Some(-start), new_kmer_set));
        }
    }
    if internal {
        let kmer_sets = kmer_chunks(adapter, max_errors + 1);
        search_sets.push((0, None, kmer_sets));
    }
    remove_redundant_kmers(&search_sets)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_chunks() {
        let result = kmer_chunks("AABCABCABC", 3);
        // Should partition into chunks of size ~3-4
        assert!(!result.is_empty());

        let result = kmer_chunks("ABCDEF", 2);
        assert!(result.contains("ABC"));
        assert!(result.contains("DEF"));

        let result = kmer_chunks("ABCDEFG", 2);
        // 7/2 = 3 remainder 1: one chunk of 4, one of 3
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_create_positions_and_kmers_back() {
        let result = create_positions_and_kmers("AGATCGGAAGAGC", 3, 0.1, true, false, false);
        assert!(!result.is_empty());
    }

    #[test]
    fn test_create_positions_and_kmers_front() {
        let result = create_positions_and_kmers("AGATCGGAAGAGC", 3, 0.1, false, true, false);
        assert!(!result.is_empty());
    }

    #[test]
    fn test_create_positions_and_kmers_both() {
        let result = create_positions_and_kmers("AGATCGGAAGAGC", 3, 0.1, true, true, true);
        assert!(!result.is_empty());
    }

    #[test]
    fn test_create_positions_and_kmers_no_errors() {
        let result = create_positions_and_kmers("ACGT", 3, 0.0, true, false, false);
        assert!(!result.is_empty());
    }
}
