/// Character encoding tables for DNA sequence matching.
///
/// Port of cutadapt/_match_tables.py.

/// Return a translation table that maps A, C, G, T characters to the lower
/// four bits of a byte. Other characters (including IUPAC) are mapped to 0x80.
///
/// Lowercase versions are also translated, and U is treated the same as T.
pub fn acgt_table() -> [u8; 256] {
    let mut t = [0x80u8; 256];
    // A=1, C=2, G=4, T=8, U=8
    t[b'A' as usize] = 1;
    t[b'a' as usize] = 1;
    t[b'C' as usize] = 2;
    t[b'c' as usize] = 2;
    t[b'G' as usize] = 4;
    t[b'g' as usize] = 4;
    t[b'T' as usize] = 8;
    t[b't' as usize] = 8;
    t[b'U' as usize] = 8;
    t[b'u' as usize] = 8;
    t
}

/// Return a translation table for IUPAC characters.
///
/// The four least significant bits represent one nucleotide each:
///   A=0x01, C=0x02, G=0x04, T=0x08
///
/// For the "N" wildcard, additionally the most significant bit is set (0x80),
/// which allows it to match characters encoded with `acgt_table`.
///
/// Whether two encoded characters x and y match can be checked with `(x & y) != 0`.
pub fn iupac_table() -> [u8; 256] {
    const A: u8 = 1;
    const C: u8 = 2;
    const G: u8 = 4;
    const T: u8 = 8;

    let mut t = [0u8; 256];

    let mappings: &[(u8, u8)] = &[
        (b'X', 0),
        (b'A', A),
        (b'C', C),
        (b'G', G),
        (b'T', T),
        (b'U', T),
        (b'R', A | G),
        (b'Y', C | T),
        (b'S', G | C),
        (b'W', A | T),
        (b'K', G | T),
        (b'M', A | C),
        (b'B', C | G | T),
        (b'D', A | G | T),
        (b'H', A | C | T),
        (b'V', A | C | G),
        (b'N', A | C | G | T | 0x80),
    ];

    for &(ch, val) in mappings {
        t[ch as usize] = val;
        t[(ch as char).to_ascii_lowercase() as usize] = val;
    }
    t
}

/// Return a table that converts ASCII to uppercase.
pub fn upper_table() -> [u8; 256] {
    let mut t = [0u8; 256];
    for i in 0..256 {
        t[i] = (i as u8 as char).to_ascii_uppercase() as u8;
    }
    t
}

/// Translate a byte sequence using a lookup table.
pub fn translate(s: &[u8], table: &[u8; 256]) -> Vec<u8> {
    s.iter().map(|&b| table[b as usize]).collect()
}

/// For each byte value in `ref_table`, generate a Vec of all byte values from
/// `query_table` that match according to `comp_op`.
///
/// `comp_op` is either equality or bitwise-AND.
fn all_matches_generator(
    ref_table: &[u8; 256],
    query_table: &[u8; 256],
    use_bitwise_and: bool,
) -> Vec<Vec<u8>> {
    (0..256)
        .map(|i| {
            let ref_char = ref_table[i];
            let mut matches = Vec::new();
            for j in 1u8..128 {
                // Only ASCII characters, skip NULL byte (j=0)
                let query_char = query_table[j as usize];
                let is_match = if use_bitwise_and {
                    (ref_char & query_char) != 0
                } else {
                    ref_char == query_char
                };
                if is_match {
                    matches.push(j);
                }
            }
            matches
        })
        .collect()
}

/// Return a lookup table for character matching.
///
/// `matches_lookup(ref_wildcards, query_wildcards)` returns a Vec of 256 entries.
/// Each entry `matches_lookup[c]` is a Vec<u8> of all characters that match `c`.
///
/// The matching mode depends on the wildcard flags:
/// - No wildcards: direct uppercase equality
/// - ref_wildcards only: IUPAC ref vs ACGT query, bitwise AND
/// - query_wildcards only: ACGT ref vs IUPAC query, bitwise AND
/// - Both wildcards: IUPAC ref vs IUPAC query, bitwise AND
pub fn matches_lookup(ref_wildcards: bool, query_wildcards: bool) -> Vec<Vec<u8>> {
    if !ref_wildcards && !query_wildcards {
        let ref_table = upper_table();
        let query_table = upper_table();
        all_matches_generator(&ref_table, &query_table, false)
    } else if ref_wildcards && !query_wildcards {
        let ref_table = iupac_table();
        let query_table = acgt_table();
        all_matches_generator(&ref_table, &query_table, true)
    } else if !ref_wildcards && query_wildcards {
        let ref_table = acgt_table();
        let query_table = iupac_table();
        all_matches_generator(&ref_table, &query_table, true)
    } else {
        let ref_table = iupac_table();
        let query_table = iupac_table();
        all_matches_generator(&ref_table, &query_table, true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_acgt_table() {
        let t = acgt_table();
        assert_eq!(t[b'A' as usize], 1);
        assert_eq!(t[b'a' as usize], 1);
        assert_eq!(t[b'C' as usize], 2);
        assert_eq!(t[b'G' as usize], 4);
        assert_eq!(t[b'T' as usize], 8);
        assert_eq!(t[b'U' as usize], 8);
        assert_eq!(t[b'N' as usize], 0x80);
        assert_eq!(t[b'X' as usize], 0x80);
    }

    #[test]
    fn test_iupac_table() {
        let t = iupac_table();
        assert_eq!(t[b'A' as usize], 0x01);
        assert_eq!(t[b'C' as usize], 0x02);
        assert_eq!(t[b'G' as usize], 0x04);
        assert_eq!(t[b'T' as usize], 0x08);
        assert_eq!(t[b'R' as usize], 0x05); // A|G
        assert_eq!(t[b'Y' as usize], 0x0A); // C|T
        assert_eq!(t[b'N' as usize], 0x8F); // A|C|G|T|0x80
        assert_eq!(t[b'X' as usize], 0x00);
        // Lowercase
        assert_eq!(t[b'a' as usize], 0x01);
        assert_eq!(t[b'n' as usize], 0x8F);
    }

    #[test]
    fn test_iupac_matching() {
        let t = iupac_table();
        // N matches everything
        assert_ne!(t[b'N' as usize] & t[b'A' as usize], 0);
        assert_ne!(t[b'N' as usize] & t[b'C' as usize], 0);
        assert_ne!(t[b'N' as usize] & t[b'G' as usize], 0);
        assert_ne!(t[b'N' as usize] & t[b'T' as usize], 0);
        // R = A|G
        assert_ne!(t[b'R' as usize] & t[b'A' as usize], 0);
        assert_ne!(t[b'R' as usize] & t[b'G' as usize], 0);
        assert_eq!(t[b'R' as usize] & t[b'C' as usize], 0);
        assert_eq!(t[b'R' as usize] & t[b'T' as usize], 0);
        // X matches nothing
        assert_eq!(t[b'X' as usize] & t[b'A' as usize], 0);
    }

    #[test]
    fn test_upper_table() {
        let t = upper_table();
        assert_eq!(t[b'a' as usize], b'A');
        assert_eq!(t[b'A' as usize], b'A');
        assert_eq!(t[b'z' as usize], b'Z');
        assert_eq!(t[b'0' as usize], b'0');
    }

    #[test]
    fn test_matches_lookup_no_wildcards() {
        let lookup = matches_lookup(false, false);
        // 'A' matches 'A' and 'a'
        assert!(lookup[b'A' as usize].contains(&b'A'));
        assert!(lookup[b'a' as usize].contains(&b'A'));
        // 'A' does not match 'C'
        assert!(!lookup[b'A' as usize].contains(&b'C'));
    }

    #[test]
    fn test_matches_lookup_ref_wildcards() {
        let lookup = matches_lookup(true, false);
        // N in reference should match A, C, G, T in query
        assert!(lookup[b'N' as usize].contains(&b'A'));
        assert!(lookup[b'N' as usize].contains(&b'C'));
        assert!(lookup[b'N' as usize].contains(&b'G'));
        assert!(lookup[b'N' as usize].contains(&b'T'));
        // R in reference should match A, G but not C, T
        assert!(lookup[b'R' as usize].contains(&b'A'));
        assert!(lookup[b'R' as usize].contains(&b'G'));
        assert!(!lookup[b'R' as usize].contains(&b'C'));
    }

    #[test]
    fn test_matches_lookup_query_wildcards() {
        let lookup = matches_lookup(false, true);
        // A in reference should match N in query (since N has A bit set)
        assert!(lookup[b'A' as usize].contains(&b'N'));
        assert!(lookup[b'A' as usize].contains(&b'R')); // R = A|G
        assert!(!lookup[b'A' as usize].contains(&b'Y')); // Y = C|T
    }

    #[test]
    fn test_matches_lookup_both_wildcards() {
        let lookup = matches_lookup(true, true);
        // N matches N
        assert!(lookup[b'N' as usize].contains(&b'N'));
        // R matches S? R=A|G, S=G|C -> overlap at G
        assert!(lookup[b'R' as usize].contains(&b'S'));
        // R matches Y? R=A|G, Y=C|T -> no overlap
        assert!(!lookup[b'R' as usize].contains(&b'Y'));
    }

    #[test]
    fn test_translate() {
        let t = acgt_table();
        let result = translate(b"ACGT", &t);
        assert_eq!(result, vec![1, 2, 4, 8]);
    }
}
