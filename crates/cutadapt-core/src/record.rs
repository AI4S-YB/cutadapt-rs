/// A sequencing read with name, sequence, and optional quality scores.
///
/// Equivalent to dnaio.SequenceRecord in the Python implementation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceRecord {
    pub name: String,
    pub sequence: String,
    pub qualities: Option<String>,
}

impl SequenceRecord {
    pub fn new(name: &str, sequence: &str, qualities: Option<&str>) -> Self {
        Self {
            name: name.to_string(),
            sequence: sequence.to_string(),
            qualities: qualities.map(|q| q.to_string()),
        }
    }

    pub fn from_parts(name: String, sequence: String, qualities: Option<String>) -> Self {
        Self {
            name,
            sequence,
            qualities,
        }
    }

    /// Return the length of the sequence.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Slice the record, returning a new record with the subsequence and sub-qualities.
    /// Equivalent to Python's `read[start:stop]`.
    pub fn slice(&self, start: usize, stop: usize) -> Self {
        let stop = stop.min(self.sequence.len());
        let start = start.min(stop);
        Self {
            name: self.name.clone(),
            sequence: self.sequence[start..stop].to_string(),
            qualities: self
                .qualities
                .as_ref()
                .map(|q| q[start..stop].to_string()),
        }
    }

    /// Return the reverse complement of this record.
    pub fn reverse_complement(&self) -> Self {
        let rc_seq: String = self
            .sequence
            .bytes()
            .rev()
            .map(complement_base)
            .map(|b| b as char)
            .collect();
        let rc_qual = self.qualities.as_ref().map(|q| {
            q.chars().rev().collect::<String>()
        });
        Self {
            name: self.name.clone(),
            sequence: rc_seq,
            qualities: rc_qual,
        }
    }
}

fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        b'N' => b'N',
        b'n' => b'n',
        other => other,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic() {
        let rec = SequenceRecord::new("read1", "ACGT", Some("IIII"));
        assert_eq!(rec.len(), 4);
        assert!(!rec.is_empty());
    }

    #[test]
    fn test_slice() {
        let rec = SequenceRecord::new("read1", "ACGTACGT", Some("IIIIIIII"));
        let sliced = rec.slice(2, 6);
        assert_eq!(sliced.sequence, "GTAC");
        assert_eq!(sliced.qualities.as_deref(), Some("IIII"));
        assert_eq!(sliced.name, "read1");
    }

    #[test]
    fn test_slice_no_qualities() {
        let rec = SequenceRecord::new("read1", "ACGTACGT", None);
        let sliced = rec.slice(2, 6);
        assert_eq!(sliced.sequence, "GTAC");
        assert_eq!(sliced.qualities, None);
    }

    #[test]
    fn test_slice_bounds() {
        let rec = SequenceRecord::new("read1", "ACGT", Some("IIII"));
        let sliced = rec.slice(2, 100);
        assert_eq!(sliced.sequence, "GT");

        let sliced = rec.slice(0, 0);
        assert_eq!(sliced.sequence, "");
    }

    #[test]
    fn test_reverse_complement() {
        let rec = SequenceRecord::new("read1", "ACGT", Some("1234"));
        let rc = rec.reverse_complement();
        assert_eq!(rc.sequence, "ACGT"); // ACGT is its own reverse complement
        assert_eq!(rc.qualities.as_deref(), Some("4321"));

        let rec2 = SequenceRecord::new("read2", "AACG", Some("1234"));
        let rc2 = rec2.reverse_complement();
        assert_eq!(rc2.sequence, "CGTT");
        assert_eq!(rc2.qualities.as_deref(), Some("4321"));
    }

    #[test]
    fn test_empty() {
        let rec = SequenceRecord::new("read1", "", None);
        assert!(rec.is_empty());
        assert_eq!(rec.len(), 0);
    }
}
