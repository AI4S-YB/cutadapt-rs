/// Read modification pipeline.
///
/// Port of cutadapt/modifiers.py.
///
/// Each modifier transforms a [`SequenceRecord`] in place, recording trimming
/// metadata in [`ModificationInfo`].

pub mod adapter;
pub mod quality;
pub mod rename;
pub mod simple;

use crate::info::ModificationInfo;
use crate::record::SequenceRecord;

// ---------------------------------------------------------------------------
// Traits
// ---------------------------------------------------------------------------

/// A modifier that operates on a single read.
pub trait SingleEndModifier {
    /// Apply the modification, returning the (possibly altered) read.
    fn modify(
        &mut self,
        read: SequenceRecord,
        info: &mut ModificationInfo,
    ) -> SequenceRecord;
}

/// A modifier that operates on a paired-end read pair.
pub trait PairedEndModifier {
    /// Apply the modification to both reads, returning the (possibly altered) pair.
    fn modify_pair(
        &mut self,
        r1: SequenceRecord,
        r2: SequenceRecord,
        info1: &mut ModificationInfo,
        info2: &mut ModificationInfo,
    ) -> (SequenceRecord, SequenceRecord);
}

// ---------------------------------------------------------------------------
// PairedEndModifierWrapper
// ---------------------------------------------------------------------------

/// Wraps up to two [`SingleEndModifier`]s so they can be applied as a
/// [`PairedEndModifier`].  At least one of the two modifiers must be present.
pub struct PairedEndModifierWrapper {
    modifier1: Option<Box<dyn SingleEndModifier>>,
    modifier2: Option<Box<dyn SingleEndModifier>>,
}

impl PairedEndModifierWrapper {
    /// Create a new wrapper.
    ///
    /// # Panics
    /// Panics if both modifiers are `None`.
    pub fn new(
        modifier1: Option<Box<dyn SingleEndModifier>>,
        modifier2: Option<Box<dyn SingleEndModifier>>,
    ) -> Self {
        assert!(
            modifier1.is_some() || modifier2.is_some(),
            "At least one modifier must be provided"
        );
        Self {
            modifier1,
            modifier2,
        }
    }
}

impl PairedEndModifier for PairedEndModifierWrapper {
    fn modify_pair(
        &mut self,
        r1: SequenceRecord,
        r2: SequenceRecord,
        info1: &mut ModificationInfo,
        info2: &mut ModificationInfo,
    ) -> (SequenceRecord, SequenceRecord) {
        let r1 = match self.modifier1 {
            Some(ref mut m) => m.modify(r1, info1),
            None => r1,
        };
        let r2 = match self.modifier2 {
            Some(ref mut m) => m.modify(r2, info2),
            None => r2,
        };
        (r1, r2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::info::ModificationInfo;
    use crate::record::SequenceRecord;

    /// A trivial modifier that uppercases the sequence (for testing).
    struct UppercaseModifier;

    impl SingleEndModifier for UppercaseModifier {
        fn modify(
            &mut self,
            mut read: SequenceRecord,
            _info: &mut ModificationInfo,
        ) -> SequenceRecord {
            read.sequence = read.sequence.to_uppercase();
            read
        }
    }

    #[test]
    fn test_paired_wrapper_both() {
        let r1 = SequenceRecord::new("r1", "acgt", Some("IIII"));
        let r2 = SequenceRecord::new("r2", "tgca", Some("IIII"));
        let mut info1 = ModificationInfo::new(r1.clone());
        let mut info2 = ModificationInfo::new(r2.clone());

        let mut wrapper = PairedEndModifierWrapper::new(
            Some(Box::new(UppercaseModifier)),
            Some(Box::new(UppercaseModifier)),
        );

        let (out1, out2) = wrapper.modify_pair(r1, r2, &mut info1, &mut info2);
        assert_eq!(out1.sequence, "ACGT");
        assert_eq!(out2.sequence, "TGCA");
    }

    #[test]
    fn test_paired_wrapper_r1_only() {
        let r1 = SequenceRecord::new("r1", "acgt", Some("IIII"));
        let r2 = SequenceRecord::new("r2", "tgca", Some("IIII"));
        let mut info1 = ModificationInfo::new(r1.clone());
        let mut info2 = ModificationInfo::new(r2.clone());

        let mut wrapper =
            PairedEndModifierWrapper::new(Some(Box::new(UppercaseModifier)), None);

        let (out1, out2) = wrapper.modify_pair(r1, r2, &mut info1, &mut info2);
        assert_eq!(out1.sequence, "ACGT");
        assert_eq!(out2.sequence, "tgca"); // untouched
    }

    #[test]
    fn test_paired_wrapper_r2_only() {
        let r1 = SequenceRecord::new("r1", "acgt", Some("IIII"));
        let r2 = SequenceRecord::new("r2", "tgca", Some("IIII"));
        let mut info1 = ModificationInfo::new(r1.clone());
        let mut info2 = ModificationInfo::new(r2.clone());

        let mut wrapper =
            PairedEndModifierWrapper::new(None, Some(Box::new(UppercaseModifier)));

        let (out1, out2) = wrapper.modify_pair(r1, r2, &mut info1, &mut info2);
        assert_eq!(out1.sequence, "acgt"); // untouched
        assert_eq!(out2.sequence, "TGCA");
    }

    #[test]
    #[should_panic(expected = "At least one modifier must be provided")]
    fn test_paired_wrapper_both_none_panics() {
        PairedEndModifierWrapper::new(None, None);
    }
}
