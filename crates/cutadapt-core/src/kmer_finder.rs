// Shift-and bitmatrix k-mer pre-screening.
//
// Port of cutadapt/_kmer_finder.pyx.

use crate::match_tables::matches_lookup;

const BITMASK_INDEX_SIZE: usize = 128;
const MAX_WORD_SIZE: usize = 64;

/// Maximum combined length of k-mers per search entry.
pub const MAXIMUM_WORD_SIZE: usize = MAX_WORD_SIZE;

/// A single search entry: one group of k-mers with shared search region.
struct KmerSearchEntry {
    mask_offset: usize,
    search_start: isize,
    search_stop: isize, // 0 means go to end of sequence
    init_mask: u64,
    found_mask: u64,
}

/// Find k-mers in strings using an enhanced shift-and algorithm.
///
/// Multiple k-mers can be packed into a single 64-bit word and searched simultaneously.
pub struct KmerFinder {
    search_entries: Vec<KmerSearchEntry>,
    search_masks: Vec<u64>, // flattened: number_of_searches * BITMASK_INDEX_SIZE
    pub ref_wildcards: bool,
    pub query_wildcards: bool,
}

impl KmerFinder {
    /// Create a new KmerFinder.
    ///
    /// `positions_and_kmers` is a list of `(start, stop, kmers)` tuples where:
    /// - `start` is the search start position (negative = from end)
    /// - `stop` is the search stop position (None/0 = end of sequence, negative = from end)
    /// - `kmers` is a list of k-mer strings to search for
    pub fn new(
        positions_and_kmers: &[(isize, Option<isize>, Vec<String>)],
        ref_wildcards: bool,
        query_wildcards: bool,
    ) -> Self {
        let match_lut = matches_lookup(ref_wildcards, query_wildcards);
        let mut search_entries = Vec::new();
        let mut search_masks: Vec<u64> = Vec::new();

        for (start, stop, kmers) in positions_and_kmers {
            let mut index = 0;
            while index < kmers.len() {
                let mut search_word = vec![0u8; 64];
                let mut offset: usize = 0;
                let mut init_mask: u64 = 0;
                let mut found_mask: u64 = 0;

                while index < kmers.len() {
                    let kmer = &kmers[index];
                    let kmer_bytes = kmer.as_bytes();
                    let kmer_length = kmer_bytes.len();
                    assert!(
                        kmer_length <= MAX_WORD_SIZE,
                        "{} of length {} is longer than the maximum of {}",
                        kmer,
                        kmer_length,
                        MAX_WORD_SIZE
                    );
                    if offset + kmer_length > MAX_WORD_SIZE {
                        break;
                    }
                    init_mask |= 1u64 << offset;
                    search_word[offset..offset + kmer_length].copy_from_slice(kmer_bytes);
                    found_mask |= 1u64 << (offset + kmer_length - 1);
                    offset += kmer_length;
                    index += 1;
                }

                let mask_offset = search_masks.len();
                // Allocate BITMASK_INDEX_SIZE entries for this search
                search_masks.resize(search_masks.len() + BITMASK_INDEX_SIZE, 0);

                // Populate needle mask
                populate_needle_mask(
                    &mut search_masks[mask_offset..mask_offset + BITMASK_INDEX_SIZE],
                    &search_word[..offset],
                    &match_lut,
                );

                let stop_val = match stop {
                    Some(s) => *s,
                    None => 0, // Encode 'end of sequence' as 0
                };

                search_entries.push(KmerSearchEntry {
                    mask_offset,
                    search_start: *start,
                    search_stop: stop_val,
                    init_mask,
                    found_mask,
                });
            }
        }

        Self {
            search_entries,
            search_masks,
            ref_wildcards,
            query_wildcards,
        }
    }

    /// Check if any of the k-mers are present in the given sequence.
    pub fn kmers_present(&self, sequence: &str) -> bool {
        let seq = sequence.as_bytes();
        let seq_length = seq.len() as isize;

        for entry in &self.search_entries {
            let mut start = entry.search_start;
            let mut stop = entry.search_stop;

            // Normalize start
            if start < 0 {
                start = seq_length + start;
                if start < 0 {
                    start = 0;
                }
            } else if start > seq_length {
                continue;
            }

            // Normalize stop
            if stop < 0 {
                stop = seq_length + stop;
                if stop <= 0 {
                    continue;
                }
            } else if stop == 0 {
                stop = seq_length;
            }

            let search_length = stop - start;
            if search_length <= 0 {
                continue;
            }

            let search_ptr = &seq[start as usize..(start + search_length) as usize];
            let mask_ptr =
                &self.search_masks[entry.mask_offset..entry.mask_offset + BITMASK_INDEX_SIZE];

            if shift_and_multiple_is_present(
                search_ptr,
                mask_ptr,
                entry.init_mask,
                entry.found_mask,
            ) {
                return true;
            }
        }
        false
    }
}

/// Populate the needle mask table for a concatenated kmer word.
fn populate_needle_mask(needle_mask: &mut [u64], needle: &[u8], match_lut: &[Vec<u8>]) {
    // needle_mask is already zeroed from resize
    for (i, &c) in needle.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let matching_chars = &match_lut[c as usize];
        for &mc in matching_chars {
            needle_mask[mc as usize] |= 1u64 << i;
        }
    }
}

/// Core shift-and algorithm for multiple pattern matching.
///
/// Returns true if any of the patterns (encoded in init_mask/found_mask)
/// are found in the haystack.
fn shift_and_multiple_is_present(
    haystack: &[u8],
    needle_mask: &[u64],
    init_mask: u64,
    found_mask: u64,
) -> bool {
    let mut r: u64 = 0;
    for &byte in haystack {
        r <<= 1;
        r |= init_mask;
        r &= needle_mask[byte as usize];
        if (r & found_mask) != 0 {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_kmer_present() {
        let kmers = vec![(0isize, None, vec!["ACGT".to_string()])];
        let finder = KmerFinder::new(&kmers, false, false);
        assert!(finder.kmers_present("XXXXXACGTXXXXX"));
        assert!(!finder.kmers_present("XXXXXACGXXXXXS"));
    }

    #[test]
    fn test_multiple_kmers() {
        let kmers = vec![(0isize, None, vec!["AAA".to_string(), "TTT".to_string()])];
        let finder = KmerFinder::new(&kmers, false, false);
        assert!(finder.kmers_present("XXXXXAAAXXX"));
        assert!(finder.kmers_present("XXXXXTTTXXX"));
        assert!(!finder.kmers_present("XXXXXCCCXXX"));
    }

    #[test]
    fn test_search_region() {
        let kmers = vec![(-5isize, None, vec!["ACG".to_string()])];
        let finder = KmerFinder::new(&kmers, false, false);
        // ACG is at position 0, but search starts at -5 from end
        assert!(!finder.kmers_present("ACGTTTTTTTT"));
        // ACG is at end, within the search region
        assert!(finder.kmers_present("TTTTTTTACGT"));
    }

    #[test]
    fn test_case_insensitive() {
        let kmers = vec![(0isize, None, vec!["ACGT".to_string()])];
        let finder = KmerFinder::new(&kmers, false, false);
        assert!(finder.kmers_present("XXXacgtXXX"));
    }

    #[test]
    fn test_wildcard_ref() {
        // N in kmer should match any base
        let kmers = vec![(0isize, None, vec!["ANGT".to_string()])];
        let finder = KmerFinder::new(&kmers, true, false);
        assert!(finder.kmers_present("XXACGTXX"));
        assert!(finder.kmers_present("XXATGTXX"));
    }

    #[test]
    fn test_wildcard_query() {
        // N in sequence should match any base in kmer
        let kmers = vec![(0isize, None, vec!["ACGT".to_string()])];
        let finder = KmerFinder::new(&kmers, false, true);
        assert!(finder.kmers_present("XXANGTXX"));
    }

    #[test]
    fn test_empty_sequence() {
        let kmers = vec![(0isize, None, vec!["ACG".to_string()])];
        let finder = KmerFinder::new(&kmers, false, false);
        assert!(!finder.kmers_present(""));
    }

    #[test]
    fn test_stop_region() {
        let kmers = vec![(0isize, Some(5), vec!["ACG".to_string()])];
        let finder = KmerFinder::new(&kmers, false, false);
        assert!(finder.kmers_present("ACGTTTTTTTT")); // ACG at pos 0, within [0,5)
        assert!(!finder.kmers_present("TTTTTACGTTT")); // ACG at pos 5, outside [0,5)
    }
}
