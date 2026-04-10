/// Simple sequence modifiers that do not depend on quality scores or adapters.

use crate::info::ModificationInfo;
use crate::record::SequenceRecord;

use super::SingleEndModifier;

// ---------------------------------------------------------------------------
// UnconditionalCutter
// ---------------------------------------------------------------------------

/// Remove a fixed number of bases from the start or end of a read.
///
/// * Positive `length` — remove from the **beginning**.
/// * Negative `length` — remove from the **end**.
///
/// Equivalent to Python's `UnconditionalCutter`.
pub struct UnconditionalCutter {
    pub length: i32,
}

impl UnconditionalCutter {
    pub fn new(length: i32) -> Self {
        Self { length }
    }
}

impl SingleEndModifier for UnconditionalCutter {
    fn modify(
        &mut self,
        read: SequenceRecord,
        info: &mut ModificationInfo,
    ) -> SequenceRecord {
        if self.length > 0 {
            let cut = (self.length as usize).min(read.len());
            info.cut_prefix = Some(read.sequence[..cut].to_string());
            read.slice(cut, read.len())
        } else if self.length < 0 {
            let abs_len = (-self.length) as usize;
            let cut_from = read.len().saturating_sub(abs_len);
            info.cut_suffix = Some(read.sequence[cut_from..].to_string());
            read.slice(0, cut_from)
        } else {
            // length == 0 → no-op
            read
        }
    }
}

// ---------------------------------------------------------------------------
// NEndTrimmer
// ---------------------------------------------------------------------------

/// Trim leading and trailing `N` characters from a read.
///
/// Equivalent to Python's `NEndTrimmer`.
pub struct NEndTrimmer;

impl NEndTrimmer {
    pub fn new() -> Self {
        Self
    }

    /// Find the first non-N position from the start.
    fn leading_n_count(seq: &[u8]) -> usize {
        seq.iter().take_while(|&&b| b == b'N' || b == b'n').count()
    }

    /// Find the first non-N position from the end.
    fn trailing_n_start(seq: &[u8]) -> usize {
        let n = seq.len();
        let trailing = seq.iter().rev().take_while(|&&b| b == b'N' || b == b'n').count();
        n - trailing
    }
}

impl Default for NEndTrimmer {
    fn default() -> Self {
        Self::new()
    }
}

impl SingleEndModifier for NEndTrimmer {
    fn modify(
        &mut self,
        read: SequenceRecord,
        _info: &mut ModificationInfo,
    ) -> SequenceRecord {
        let seq = read.sequence.as_bytes();
        let start = Self::leading_n_count(seq);
        let stop = Self::trailing_n_start(seq);
        if start >= stop {
            read.slice(0, 0)
        } else {
            read.slice(start, stop)
        }
    }
}

// ---------------------------------------------------------------------------
// Shortener
// ---------------------------------------------------------------------------

/// Unconditionally truncate a read to at most `length` bases.
///
/// * Positive `length` — keep the first `length` bases (trim from end).
/// * Negative `length` — keep the last `|length|` bases (trim from start).
///
/// Equivalent to Python's `Shortener`.
pub struct Shortener {
    pub length: i32,
}

impl Shortener {
    pub fn new(length: i32) -> Self {
        Self { length }
    }
}

impl SingleEndModifier for Shortener {
    fn modify(
        &mut self,
        read: SequenceRecord,
        _info: &mut ModificationInfo,
    ) -> SequenceRecord {
        if self.length >= 0 {
            let stop = (self.length as usize).min(read.len());
            read.slice(0, stop)
        } else {
            let abs_len = (-self.length) as usize;
            let start = read.len().saturating_sub(abs_len);
            read.slice(start, read.len())
        }
    }
}

// ---------------------------------------------------------------------------
// ZeroCapper
// ---------------------------------------------------------------------------

/// Replace negative quality values (below the encoding base) with zero.
///
/// Equivalent to Python's `ZeroCapper`.
pub struct ZeroCapper {
    pub quality_base: u8,
}

impl ZeroCapper {
    pub fn new(quality_base: u8) -> Self {
        Self { quality_base }
    }
}

impl SingleEndModifier for ZeroCapper {
    fn modify(
        &mut self,
        mut read: SequenceRecord,
        _info: &mut ModificationInfo,
    ) -> SequenceRecord {
        if let Some(ref quals) = read.qualities {
            let base = self.quality_base;
            let capped: String = quals
                .bytes()
                .map(|b| if b < base { base as char } else { b as char })
                .collect();
            read.qualities = Some(capped);
        }
        read
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::info::ModificationInfo;
    use crate::record::SequenceRecord;

    // -- UnconditionalCutter --

    #[test]
    fn test_unconditional_cutter_positive() {
        let read = SequenceRecord::new("r1", "ACGTACGT", Some("IIIIIIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut cutter = UnconditionalCutter::new(3);

        let result = cutter.modify(read, &mut info);
        assert_eq!(result.sequence, "TACGT");
        assert_eq!(info.cut_prefix.as_deref(), Some("ACG"));
    }

    #[test]
    fn test_unconditional_cutter_negative() {
        let read = SequenceRecord::new("r1", "ACGTACGT", Some("IIIIIIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut cutter = UnconditionalCutter::new(-3);

        let result = cutter.modify(read, &mut info);
        assert_eq!(result.sequence, "ACGTA");
        assert_eq!(info.cut_suffix.as_deref(), Some("CGT"));
    }

    #[test]
    fn test_unconditional_cutter_zero() {
        let read = SequenceRecord::new("r1", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut cutter = UnconditionalCutter::new(0);

        let result = cutter.modify(read.clone(), &mut info);
        assert_eq!(result.sequence, read.sequence);
    }

    #[test]
    fn test_unconditional_cutter_longer_than_read() {
        let read = SequenceRecord::new("r1", "AC", Some("II"));
        let mut info = ModificationInfo::new(read.clone());
        let mut cutter = UnconditionalCutter::new(10);

        let result = cutter.modify(read, &mut info);
        assert_eq!(result.sequence, "");
        assert_eq!(info.cut_prefix.as_deref(), Some("AC"));
    }

    // -- NEndTrimmer --

    #[test]
    fn test_n_end_trimmer_both_ends() {
        let read = SequenceRecord::new("r1", "NNNACGTNNN", Some("IIIIIIIIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut trimmer = NEndTrimmer::new();

        let result = trimmer.modify(read, &mut info);
        assert_eq!(result.sequence, "ACGT");
    }

    #[test]
    fn test_n_end_trimmer_no_ns() {
        let read = SequenceRecord::new("r1", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut trimmer = NEndTrimmer::new();

        let result = trimmer.modify(read, &mut info);
        assert_eq!(result.sequence, "ACGT");
    }

    #[test]
    fn test_n_end_trimmer_all_ns() {
        let read = SequenceRecord::new("r1", "NNNN", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut trimmer = NEndTrimmer::new();

        let result = trimmer.modify(read, &mut info);
        assert_eq!(result.sequence, "");
    }

    #[test]
    fn test_n_end_trimmer_leading_only() {
        let read = SequenceRecord::new("r1", "NNACGT", Some("IIIIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut trimmer = NEndTrimmer::new();

        let result = trimmer.modify(read, &mut info);
        assert_eq!(result.sequence, "ACGT");
    }

    #[test]
    fn test_n_end_trimmer_trailing_only() {
        let read = SequenceRecord::new("r1", "ACGTNN", Some("IIIIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut trimmer = NEndTrimmer::new();

        let result = trimmer.modify(read, &mut info);
        assert_eq!(result.sequence, "ACGT");
    }

    // -- Shortener --

    #[test]
    fn test_shortener_positive() {
        let read = SequenceRecord::new("r1", "ACGTACGT", Some("IIIIIIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut shortener = Shortener::new(4);

        let result = shortener.modify(read, &mut info);
        assert_eq!(result.sequence, "ACGT");
    }

    #[test]
    fn test_shortener_longer_than_read() {
        let read = SequenceRecord::new("r1", "AC", Some("II"));
        let mut info = ModificationInfo::new(read.clone());
        let mut shortener = Shortener::new(10);

        let result = shortener.modify(read, &mut info);
        assert_eq!(result.sequence, "AC");
    }

    #[test]
    fn test_shortener_negative() {
        let read = SequenceRecord::new("r1", "ACGTACGT", Some("IIIIIIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut shortener = Shortener::new(-4);

        let result = shortener.modify(read, &mut info);
        assert_eq!(result.sequence, "ACGT");
        assert_eq!(result.qualities.as_deref(), Some("IIII"));
    }

    #[test]
    fn test_shortener_zero() {
        let read = SequenceRecord::new("r1", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut shortener = Shortener::new(0);

        let result = shortener.modify(read, &mut info);
        assert_eq!(result.sequence, "");
    }

    // -- ZeroCapper --

    #[test]
    fn test_zero_capper() {
        // Quality characters below base 33 should be capped to base
        let read = SequenceRecord::new("r1", "ACGT", Some("\x20!I\x1F"));
        // \x20 = 32 < 33 -> should become 33 ('!')
        // '!' = 33 -> stays
        // 'I' = 73 -> stays
        // \x1F = 31 < 33 -> should become 33 ('!')
        let mut info = ModificationInfo::new(read.clone());
        let mut capper = ZeroCapper::new(33);

        let result = capper.modify(read, &mut info);
        assert_eq!(result.qualities.as_deref(), Some("!!I!"));
    }

    #[test]
    fn test_zero_capper_no_qualities() {
        let read = SequenceRecord::new("r1", "ACGT", None);
        let mut info = ModificationInfo::new(read.clone());
        let mut capper = ZeroCapper::new(33);

        let result = capper.modify(read.clone(), &mut info);
        assert_eq!(result.qualities, None);
    }

    #[test]
    fn test_zero_capper_all_good() {
        let read = SequenceRecord::new("r1", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut capper = ZeroCapper::new(33);

        let result = capper.modify(read, &mut info);
        assert_eq!(result.qualities.as_deref(), Some("IIII"));
    }
}
