/// Quality-based trimming modifiers.
///
/// Wraps the low-level functions in [`crate::qualtrim`].

use std::collections::HashMap;

use crate::info::ModificationInfo;
use crate::qualtrim;
use crate::record::SequenceRecord;

use super::SingleEndModifier;

// ---------------------------------------------------------------------------
// QualityTrimmer
// ---------------------------------------------------------------------------

/// Trim low-quality bases from both ends of a read using the BWA algorithm.
///
/// Equivalent to Python's `QualityTrimmer`.
pub struct QualityTrimmer {
    pub cutoff_front: i32,
    pub cutoff_back: i32,
    pub base: i32,
    /// Running total of bases removed.
    pub trimmed_bases: usize,
}

impl QualityTrimmer {
    pub fn new(cutoff_front: i32, cutoff_back: i32, base: i32) -> Self {
        Self {
            cutoff_front,
            cutoff_back,
            base,
            trimmed_bases: 0,
        }
    }
}

impl SingleEndModifier for QualityTrimmer {
    fn modify(
        &mut self,
        read: SequenceRecord,
        _info: &mut ModificationInfo,
    ) -> SequenceRecord {
        let qualities = match read.qualities.as_ref() {
            Some(q) => q.as_bytes(),
            None => return read,
        };
        let (start, stop) =
            qualtrim::quality_trim_index(qualities, self.cutoff_front, self.cutoff_back, self.base);
        self.trimmed_bases += read.len() - (stop - start);
        read.slice(start, stop)
    }
}

// ---------------------------------------------------------------------------
// NextseqQualityTrimmer
// ---------------------------------------------------------------------------

/// NextSeq-specific quality trimmer that also trims trailing high-quality G
/// bases (dark cycles).
///
/// Equivalent to Python's `NextseqQualityTrimmer`.
pub struct NextseqQualityTrimmer {
    pub cutoff: i32,
    pub base: i32,
    /// Running total of bases removed.
    pub trimmed_bases: usize,
}

impl NextseqQualityTrimmer {
    pub fn new(cutoff: i32, base: i32) -> Self {
        Self {
            cutoff,
            base,
            trimmed_bases: 0,
        }
    }
}

impl SingleEndModifier for NextseqQualityTrimmer {
    fn modify(
        &mut self,
        read: SequenceRecord,
        _info: &mut ModificationInfo,
    ) -> SequenceRecord {
        let qualities = match read.qualities.as_ref() {
            Some(q) => q.as_bytes(),
            None => return read,
        };
        let stop = qualtrim::nextseq_trim_index(
            read.sequence.as_bytes(),
            qualities,
            self.cutoff,
            self.base,
        );
        self.trimmed_bases += read.len() - stop;
        read.slice(0, stop)
    }
}

// ---------------------------------------------------------------------------
// PolyATrimmer
// ---------------------------------------------------------------------------

/// Trim poly-A tails (or poly-T heads when `revcomp` is set).
///
/// Equivalent to Python's `PolyATrimmer`.
pub struct PolyATrimmer {
    pub revcomp: bool,
    /// Histogram: number of bases trimmed -> count.
    pub trimmed_bases: HashMap<usize, usize>,
}

impl PolyATrimmer {
    pub fn new(revcomp: bool) -> Self {
        Self {
            revcomp,
            trimmed_bases: HashMap::new(),
        }
    }
}

impl SingleEndModifier for PolyATrimmer {
    fn modify(
        &mut self,
        read: SequenceRecord,
        _info: &mut ModificationInfo,
    ) -> SequenceRecord {
        if self.revcomp {
            let index = qualtrim::poly_a_trim_index(read.sequence.as_bytes(), true);
            *self.trimmed_bases.entry(index).or_insert(0) += 1;
            read.slice(index, read.len())
        } else {
            let index = qualtrim::poly_a_trim_index(read.sequence.as_bytes(), false);
            let trimmed = read.len() - index;
            *self.trimmed_bases.entry(trimmed).or_insert(0) += 1;
            read.slice(0, index)
        }
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

    // -- QualityTrimmer --

    #[test]
    fn test_quality_trimmer_trims_both_ends() {
        let read = SequenceRecord::new("r1", "AACCGGTT", Some("##IIII##"));
        let mut info = ModificationInfo::new(read.clone());
        let mut trimmer = QualityTrimmer::new(10, 10, 33);

        let result = trimmer.modify(read, &mut info);
        assert_eq!(result.sequence, "CCGG");
        assert_eq!(result.qualities.as_deref(), Some("IIII"));
        assert_eq!(trimmer.trimmed_bases, 4);
    }

    #[test]
    fn test_quality_trimmer_no_qualities() {
        let read = SequenceRecord::new("r1", "ACGT", None);
        let mut info = ModificationInfo::new(read.clone());
        let mut trimmer = QualityTrimmer::new(10, 10, 33);

        let result = trimmer.modify(read.clone(), &mut info);
        assert_eq!(result.sequence, read.sequence);
        assert_eq!(trimmer.trimmed_bases, 0);
    }

    #[test]
    fn test_quality_trimmer_no_trim_needed() {
        let read = SequenceRecord::new("r1", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut trimmer = QualityTrimmer::new(10, 10, 33);

        let result = trimmer.modify(read, &mut info);
        assert_eq!(result.sequence, "ACGT");
        assert_eq!(trimmer.trimmed_bases, 0);
    }

    #[test]
    fn test_quality_trimmer_accumulates() {
        let mut trimmer = QualityTrimmer::new(10, 10, 33);

        let read1 = SequenceRecord::new("r1", "AACCGGTT", Some("##IIII##"));
        let mut info1 = ModificationInfo::new(read1.clone());
        trimmer.modify(read1, &mut info1);

        let read2 = SequenceRecord::new("r2", "AACCGGTT", Some("##IIII##"));
        let mut info2 = ModificationInfo::new(read2.clone());
        trimmer.modify(read2, &mut info2);

        assert_eq!(trimmer.trimmed_bases, 8);
    }

    // -- NextseqQualityTrimmer --

    #[test]
    fn test_nextseq_trimmer_trims_trailing_g() {
        let read = SequenceRecord::new("r1", "ACGTGGG", Some("IIIIIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut trimmer = NextseqQualityTrimmer::new(10, 33);

        let result = trimmer.modify(read, &mut info);
        assert_eq!(result.sequence, "ACGT");
        assert_eq!(trimmer.trimmed_bases, 3);
    }

    #[test]
    fn test_nextseq_trimmer_no_g() {
        let read = SequenceRecord::new("r1", "ACGTACT", Some("IIIIIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut trimmer = NextseqQualityTrimmer::new(10, 33);

        let result = trimmer.modify(read, &mut info);
        assert_eq!(result.sequence, "ACGTACT");
        assert_eq!(trimmer.trimmed_bases, 0);
    }

    #[test]
    fn test_nextseq_trimmer_no_qualities() {
        let read = SequenceRecord::new("r1", "ACGTGGG", None);
        let mut info = ModificationInfo::new(read.clone());
        let mut trimmer = NextseqQualityTrimmer::new(10, 33);

        let result = trimmer.modify(read, &mut info);
        assert_eq!(result.sequence, "ACGTGGG");
        assert_eq!(trimmer.trimmed_bases, 0);
    }

    // -- PolyATrimmer --

    #[test]
    fn test_poly_a_trimmer_trims_tail() {
        let read = SequenceRecord::new("r1", "ACGTAAAA", Some("IIIIIIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut trimmer = PolyATrimmer::new(false);

        let result = trimmer.modify(read, &mut info);
        assert_eq!(result.sequence, "ACGT");
        assert_eq!(trimmer.trimmed_bases.get(&4), Some(&1));
    }

    #[test]
    fn test_poly_a_trimmer_short_tail_ignored() {
        let read = SequenceRecord::new("r1", "ACGTAA", Some("IIIIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut trimmer = PolyATrimmer::new(false);

        let result = trimmer.modify(read, &mut info);
        assert_eq!(result.sequence, "ACGTAA");
        assert_eq!(trimmer.trimmed_bases.get(&0), Some(&1));
    }

    #[test]
    fn test_poly_t_trimmer_revcomp() {
        let read = SequenceRecord::new("r1", "TTTTACGT", Some("IIIIIIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut trimmer = PolyATrimmer::new(true);

        let result = trimmer.modify(read, &mut info);
        assert_eq!(result.sequence, "ACGT");
        assert_eq!(trimmer.trimmed_bases.get(&4), Some(&1));
    }
}
