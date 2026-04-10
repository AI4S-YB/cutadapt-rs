/// Quality trimming algorithms.
///
/// Port of cutadapt/qualtrim.pyx.

use thiserror::Error;

#[derive(Debug, Error)]
pub enum QualTrimError {
    #[error("Cannot do quality trimming when no qualities are available")]
    NoQualities,
    #[error("Not a valid phred value {value} for character {character}")]
    InvalidPhred { value: u8, character: char },
}

/// Find the positions at which to trim low-quality ends from a nucleotide sequence.
/// Returns `(start, stop)` indicating the good-quality segment.
///
/// Qualities are assumed to be ASCII-encoded as `chr(qual + base)`.
///
/// The algorithm is the same as the one used by BWA within the function
/// `bwa_trim_read`:
/// - Subtract the cutoff value from all qualities.
/// - Compute partial sums from all indices to the end of the sequence.
/// - Trim sequence at the index at which the sum is minimal.
pub fn quality_trim_index(
    qualities: &[u8],
    cutoff_front: i32,
    cutoff_back: i32,
    base: i32,
) -> (usize, usize) {
    let n = qualities.len();
    let mut stop = n;
    let mut start = 0usize;

    // Find trim position for 5' end
    let mut s: i32 = 0;
    let mut max_qual: i32 = 0;
    for i in 0..n {
        s += cutoff_front - (qualities[i] as i32 - base);
        if s < 0 {
            break;
        }
        if s > max_qual {
            max_qual = s;
            start = i + 1;
        }
    }

    // Same for 3' end
    max_qual = 0;
    s = 0;
    for i in (0..n).rev() {
        s += cutoff_back - (qualities[i] as i32 - base);
        if s < 0 {
            break;
        }
        if s > max_qual {
            max_qual = s;
            stop = i;
        }
    }

    if start >= stop {
        (0, 0)
    } else {
        (start, stop)
    }
}

/// Variant of quality trimming for NextSeq data.
///
/// With Illumina NextSeq, bases are encoded with two colors. 'No color' (a dark cycle)
/// usually means that a 'G' was sequenced, but that also occurs when sequencing falls
/// off the end of the fragment. The read then contains a run of high-quality G bases
/// in the end.
///
/// This routine works as the one above, but counts qualities belonging to 'G'
/// bases as being equal to `cutoff - 1`.
pub fn nextseq_trim_index(
    sequence: &[u8],
    qualities: &[u8],
    cutoff: i32,
    base: i32,
) -> usize {
    let n = qualities.len();
    let mut s: i32 = 0;
    let mut max_qual: i32 = 0;
    let mut max_i = n;

    for i in (0..n).rev() {
        let mut q = qualities[i] as i32 - base;
        if sequence[i] == b'G' {
            q = cutoff - 1;
        }
        s += cutoff - q;
        if s < 0 {
            break;
        }
        if s > max_qual {
            max_qual = s;
            max_i = i;
        }
    }
    max_i
}

/// Return start index of poly-A tail.
///
/// If `revcomp` is true, return end of poly-T head instead.
///
/// Poly-A tails shorter than 3 are ignored.
pub fn poly_a_trim_index(sequence: &[u8], revcomp: bool) -> usize {
    let n = sequence.len();
    let mut best_score: i32 = 0;
    let mut score: i32 = 0;
    let mut errors: i32 = 0;

    if revcomp {
        // Scan forward for poly-T head
        let mut best_index: usize = 0;
        for i in 0..n {
            if sequence[i] == b'T' {
                score += 1;
            } else {
                score -= 2;
                errors += 1;
            }
            if score > best_score && errors * 5 <= (i as i32 + 1) {
                best_score = score;
                best_index = i + 1;
            }
        }
        if best_index < 3 {
            best_index = 0;
        }
        best_index
    } else {
        // Scan backward for poly-A tail
        let mut best_index: usize = n;
        for i in (0..n).rev() {
            if sequence[i] == b'A' {
                score += 1;
            } else {
                score -= 2;
                errors += 1;
            }
            if score > best_score && errors * 5 <= (n as i32 - i as i32) {
                best_score = score;
                best_index = i;
            }
        }
        if best_index > n - 3 {
            best_index = n;
        }
        best_index
    }
}

/// Lookup table: Phred score -> error probability.
///
/// `SCORE_TO_ERROR_RATE[q] = 10^(-q/10)`
const SCORE_TO_ERROR_RATE: [f64; 94] = [
    1.0,                     // 0
    0.7943282347242815,      // 1
    0.6309573444801932,      // 2
    0.5011872336272722,      // 3
    0.3981071705534972,      // 4
    0.31622776601683794,     // 5
    0.251188643150958,       // 6
    0.19952623149688797,     // 7
    0.15848931924611134,     // 8
    0.12589254117941673,     // 9
    0.1,                     // 10
    0.07943282347242814,     // 11
    0.06309573444801933,     // 12
    0.05011872336272722,     // 13
    0.039810717055349734,    // 14
    0.03162277660168379,     // 15
    0.025118864315095794,    // 16
    0.0199526231496888,      // 17
    0.015848931924611134,    // 18
    0.012589254117941675,    // 19
    0.01,                    // 20
    0.007943282347242814,    // 21
    0.00630957344480193,     // 22
    0.005011872336272725,    // 23
    0.003981071705534973,    // 24
    0.0031622776601683794,   // 25
    0.0025118864315095794,   // 26
    0.001995262314968879,    // 27
    0.001584893192461114,    // 28
    0.0012589254117941675,   // 29
    0.001,                   // 30
    0.0007943282347242813,   // 31
    0.000630957344480193,    // 32
    0.0005011872336272725,   // 33
    0.00039810717055349735,  // 34
    0.00031622776601683794,  // 35
    0.00025118864315095795,  // 36
    0.00019952623149688788,  // 37
    0.00015848931924611142,  // 38
    0.00012589254117941674,  // 39
    0.0001,                  // 40
    7.943282347242822e-05,   // 41
    6.309573444801929e-05,   // 42
    5.011872336272725e-05,   // 43
    3.9810717055349695e-05,  // 44
    3.1622776601683795e-05,  // 45
    2.5118864315095822e-05,  // 46
    1.9952623149688786e-05,  // 47
    1.584893192461114e-05,   // 48
    1.2589254117941661e-05,  // 49
    1e-05,                   // 50
    7.943282347242822e-06,   // 51
    6.30957344480193e-06,    // 52
    5.011872336272725e-06,   // 53
    3.981071705534969e-06,   // 54
    3.162277660168379e-06,   // 55
    2.5118864315095823e-06,  // 56
    1.9952623149688787e-06,  // 57
    1.584893192461114e-06,   // 58
    1.2589254117941661e-06,  // 59
    1e-06,                   // 60
    7.943282347242822e-07,   // 61
    6.30957344480193e-07,    // 62
    5.011872336272725e-07,   // 63
    3.981071705534969e-07,   // 64
    3.162277660168379e-07,   // 65
    2.5118864315095823e-07,  // 66
    1.9952623149688787e-07,  // 67
    1.584893192461114e-07,   // 68
    1.2589254117941662e-07,  // 69
    1e-07,                   // 70
    7.943282347242822e-08,   // 71
    6.30957344480193e-08,    // 72
    5.011872336272725e-08,   // 73
    3.981071705534969e-08,   // 74
    3.162277660168379e-08,   // 75
    2.511886431509582e-08,   // 76
    1.9952623149688786e-08,  // 77
    1.5848931924611143e-08,  // 78
    1.2589254117941661e-08,  // 79
    1e-08,                   // 80
    7.943282347242822e-09,   // 81
    6.309573444801943e-09,   // 82
    5.011872336272715e-09,   // 83
    3.981071705534969e-09,   // 84
    3.1622776601683795e-09,  // 85
    2.511886431509582e-09,   // 86
    1.9952623149688828e-09,  // 87
    1.584893192461111e-09,   // 88
    1.2589254117941663e-09,  // 89
    1e-09,                   // 90
    7.943282347242822e-10,   // 91
    6.309573444801942e-10,   // 92
    5.011872336272714e-10,   // 93
];

/// Return the number of expected errors from a read's qualities.
///
/// Uses the formula in Edgar et al. (2015): `sum(10^(-q/10))` for each quality score.
///
/// `qualities` are ASCII-encoded as `chr(qual + base)`.
pub fn expected_errors(qualities: &[u8], base: u8) -> Result<f64, QualTrimError> {
    let max_phred = 126u8.saturating_sub(base);
    let n = qualities.len();

    // Unrolled loop matching the C implementation
    let mut ee0: f64 = 0.0;
    let mut ee1: f64 = 0.0;
    let mut ee2: f64 = 0.0;
    let mut ee3: f64 = 0.0;

    let unroll_end = if n >= 3 { n - 3 } else { 0 };
    let mut i = 0;

    while i < unroll_end {
        let p0 = qualities[i].wrapping_sub(base);
        let p1 = qualities[i + 1].wrapping_sub(base);
        let p2 = qualities[i + 2].wrapping_sub(base);
        let p3 = qualities[i + 3].wrapping_sub(base);
        if p0 > max_phred || p1 > max_phred || p2 > max_phred || p3 > max_phred {
            // Find the exact invalid character
            for &q in qualities {
                let phred = q.wrapping_sub(base);
                if phred > max_phred {
                    return Err(QualTrimError::InvalidPhred {
                        value: q,
                        character: q as char,
                    });
                }
            }
        }
        ee0 += SCORE_TO_ERROR_RATE[p0 as usize];
        ee1 += SCORE_TO_ERROR_RATE[p1 as usize];
        ee2 += SCORE_TO_ERROR_RATE[p2 as usize];
        ee3 += SCORE_TO_ERROR_RATE[p3 as usize];
        i += 4;
    }

    while i < n {
        let phred = qualities[i].wrapping_sub(base);
        if phred > max_phred {
            return Err(QualTrimError::InvalidPhred {
                value: qualities[i],
                character: qualities[i] as char,
            });
        }
        ee0 += SCORE_TO_ERROR_RATE[phred as usize];
        i += 1;
    }

    Ok(ee0 + ee1 + ee2 + ee3)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quality_trim_both_ends() {
        // Quality string where ends are low quality
        let quals = b"##IIIIII##";
        let (start, stop) = quality_trim_index(quals, 10, 10, 33);
        assert_eq!(start, 2);
        assert_eq!(stop, 8);
    }

    #[test]
    fn test_quality_trim_no_trim_needed() {
        let quals = b"IIIIII";
        let (start, stop) = quality_trim_index(quals, 10, 10, 33);
        assert_eq!(start, 0);
        assert_eq!(stop, 6);
    }

    #[test]
    fn test_quality_trim_all_bad() {
        let quals = b"####";
        let (start, stop) = quality_trim_index(quals, 30, 30, 33);
        assert_eq!(start, 0);
        assert_eq!(stop, 0);
    }

    #[test]
    fn test_quality_trim_empty() {
        let quals = b"";
        let (start, stop) = quality_trim_index(quals, 10, 10, 33);
        assert_eq!(start, 0);
        assert_eq!(stop, 0);
    }

    #[test]
    fn test_quality_trim_front_only() {
        let quals = b"##IIIIII";
        let (start, stop) = quality_trim_index(quals, 10, 0, 33);
        assert_eq!(start, 2);
        assert_eq!(stop, 8);
    }

    #[test]
    fn test_quality_trim_back_only() {
        let quals = b"IIIIII##";
        let (start, stop) = quality_trim_index(quals, 0, 10, 33);
        assert_eq!(start, 0);
        assert_eq!(stop, 6);
    }

    #[test]
    fn test_nextseq_trim() {
        // G bases at end with high quality should still be trimmed
        let seq = b"ACGTGGG";
        let qual = b"IIIIIII"; // All high quality
        let idx = nextseq_trim_index(seq, qual, 10, 33);
        assert_eq!(idx, 4); // Trim off the three G's
    }

    #[test]
    fn test_nextseq_trim_no_g() {
        let seq = b"ACGTACT";
        let qual = b"IIIIIII";
        let idx = nextseq_trim_index(seq, qual, 10, 33);
        assert_eq!(idx, 7); // No trimming
    }

    #[test]
    fn test_poly_a_trim() {
        let seq = b"ACGTAAAA";
        let idx = poly_a_trim_index(seq, false);
        assert_eq!(idx, 4);
    }

    #[test]
    fn test_poly_a_trim_short() {
        // Poly-A tails shorter than 3 are ignored
        let seq = b"ACGTAA";
        let idx = poly_a_trim_index(seq, false);
        assert_eq!(idx, 6); // No trimming (only 2 A's)
    }

    #[test]
    fn test_poly_t_trim_revcomp() {
        let seq = b"TTTTACGT";
        let idx = poly_a_trim_index(seq, true);
        assert_eq!(idx, 4);
    }

    #[test]
    fn test_poly_t_trim_short_revcomp() {
        let seq = b"TTACGT";
        let idx = poly_a_trim_index(seq, true);
        assert_eq!(idx, 0); // Only 2 T's, below minimum 3
    }

    #[test]
    fn test_expected_errors() {
        // All quality 40 (phred 40 = very low error rate)
        let quals = b"IIIII"; // 'I' = 73, phred = 73-33 = 40
        let ee = expected_errors(quals, 33).unwrap();
        assert!(ee < 0.001);

        // All quality 0 (phred 0 = error rate 1.0)
        let quals = b"!!!!!"; // '!' = 33, phred = 33-33 = 0
        let ee = expected_errors(quals, 33).unwrap();
        assert!((ee - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_expected_errors_empty() {
        let quals = b"";
        let ee = expected_errors(quals, 33).unwrap();
        assert_eq!(ee, 0.0);
    }

    #[test]
    fn test_expected_errors_single() {
        // Phred 10 -> error rate 0.1
        let quals = &[33 + 10]; // phred 10
        let ee = expected_errors(quals, 33).unwrap();
        assert!((ee - 0.1).abs() < 1e-10);
    }
}
