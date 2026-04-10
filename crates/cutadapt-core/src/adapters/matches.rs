// Match types for adapter trimming results.
//
// Port of the Match / SingleMatch / RemoveBeforeMatch / RemoveAfterMatch / LinkedMatch
// classes from cutadapt/adapters.py.

use crate::record::SequenceRecord;

/// Common fields for a single adapter match.
#[derive(Debug, Clone)]
pub struct SingleMatchFields {
    /// Start position in adapter (reference).
    pub astart: usize,
    /// Stop position in adapter (reference).
    pub astop: usize,
    /// Start position in read (query).
    pub rstart: usize,
    /// Stop position in read (query).
    pub rstop: usize,
    /// Alignment score.
    pub score: i32,
    /// Number of errors in the alignment.
    pub errors: i32,
    /// The read sequence that was matched against.
    pub sequence: String,
    /// Number of aligned characters in the adapter. May differ from
    /// aligned characters in the read when indels are present.
    pub length: usize,
    /// Name of the adapter that produced this match.
    pub adapter_name: String,
    /// The adapter sequence itself.
    pub adapter_sequence: String,
}

impl SingleMatchFields {
    pub fn new(
        astart: usize,
        astop: usize,
        rstart: usize,
        rstop: usize,
        score: i32,
        errors: i32,
        sequence: &str,
        adapter_name: &str,
        adapter_sequence: &str,
    ) -> Self {
        Self {
            astart,
            astop,
            rstart,
            rstop,
            score,
            errors,
            sequence: sequence.to_string(),
            length: astop - astart,
            adapter_name: adapter_name.to_string(),
            adapter_sequence: adapter_sequence.to_string(),
        }
    }

    /// Return the wildcard characters matched in this alignment.
    pub fn wildcards(&self, wildcard_char: u8) -> String {
        let seq_bytes = self.sequence.as_bytes();
        let adapter_bytes = self.adapter_sequence.as_bytes();
        let mut result = Vec::new();
        for i in 0..self.length {
            let adapter_idx = self.astart + i;
            let seq_idx = self.rstart + i;
            if adapter_idx < adapter_bytes.len()
                && seq_idx < seq_bytes.len()
                && adapter_bytes[adapter_idx] == wildcard_char
            {
                result.push(seq_bytes[seq_idx]);
            }
        }
        String::from_utf8(result).unwrap_or_default()
    }

    /// Return the portion of the read that aligned to the adapter.
    pub fn match_sequence(&self) -> &str {
        &self.sequence[self.rstart..self.rstop]
    }

    /// Build an info record for this match.
    pub fn get_info_record(&self, read: &SequenceRecord) -> Vec<String> {
        let seq = &read.sequence;
        let quals = read.qualities.as_deref().unwrap_or("");
        let mut info = vec![
            String::new(),
            format!("{}", self.errors),
            format!("{}", self.rstart),
            format!("{}", self.rstop),
            seq[..self.rstart].to_string(),
            seq[self.rstart..self.rstop].to_string(),
            seq[self.rstop..].to_string(),
            self.adapter_name.clone(),
        ];
        if !quals.is_empty() {
            info.push(quals[..self.rstart].to_string());
            info.push(quals[self.rstart..self.rstop].to_string());
            info.push(quals[self.rstop..].to_string());
        } else {
            info.push(String::new());
            info.push(String::new());
            info.push(String::new());
        }
        info
    }
}

/// A match where sequence before the match is removed (5' / front adapters).
#[derive(Debug, Clone)]
pub struct RemoveBeforeMatch {
    pub inner: SingleMatchFields,
}

impl RemoveBeforeMatch {
    pub fn new(
        astart: usize,
        astop: usize,
        rstart: usize,
        rstop: usize,
        score: i32,
        errors: i32,
        sequence: &str,
        adapter_name: &str,
        adapter_sequence: &str,
    ) -> Self {
        Self {
            inner: SingleMatchFields::new(
                astart,
                astop,
                rstart,
                rstop,
                score,
                errors,
                sequence,
                adapter_name,
                adapter_sequence,
            ),
        }
    }

    /// The part of the read before the match (the removed prefix).
    pub fn rest(&self) -> &str {
        &self.inner.sequence[..self.inner.rstart]
    }

    /// Interval of the read that remains after trimming: `[rstop, len)`.
    pub fn remainder_interval(&self) -> (usize, usize) {
        (self.inner.rstop, self.inner.sequence.len())
    }

    /// Interval including retained adapter portion.
    pub fn retained_adapter_interval(&self) -> (usize, usize) {
        (self.inner.rstart, self.inner.sequence.len())
    }

    /// Trim the read, keeping only the part after the match.
    pub fn trimmed(&self, read: &SequenceRecord) -> SequenceRecord {
        read.slice(self.inner.rstop, read.len())
    }

    /// Length of the removed sequence.
    pub fn removed_sequence_length(&self) -> usize {
        self.inner.rstop
    }

    /// Convenience accessors
    pub fn score(&self) -> i32 {
        self.inner.score
    }
    pub fn errors(&self) -> i32 {
        self.inner.errors
    }
    pub fn rstart(&self) -> usize {
        self.inner.rstart
    }
    pub fn rstop(&self) -> usize {
        self.inner.rstop
    }

    pub fn match_sequence(&self) -> &str {
        self.inner.match_sequence()
    }

    pub fn get_info_records(&self, read: &SequenceRecord) -> Vec<Vec<String>> {
        vec![self.inner.get_info_record(read)]
    }
}

/// A match where sequence after the match is removed (3' / back adapters).
#[derive(Debug, Clone)]
pub struct RemoveAfterMatch {
    pub inner: SingleMatchFields,
}

impl RemoveAfterMatch {
    pub fn new(
        astart: usize,
        astop: usize,
        rstart: usize,
        rstop: usize,
        score: i32,
        errors: i32,
        sequence: &str,
        adapter_name: &str,
        adapter_sequence: &str,
    ) -> Self {
        Self {
            inner: SingleMatchFields::new(
                astart,
                astop,
                rstart,
                rstop,
                score,
                errors,
                sequence,
                adapter_name,
                adapter_sequence,
            ),
        }
    }

    /// The part of the read after the match (the removed suffix).
    pub fn rest(&self) -> &str {
        &self.inner.sequence[self.inner.rstop..]
    }

    /// Interval of the read that remains after trimming: `[0, rstart)`.
    pub fn remainder_interval(&self) -> (usize, usize) {
        (0, self.inner.rstart)
    }

    /// Interval including retained adapter portion.
    pub fn retained_adapter_interval(&self) -> (usize, usize) {
        (0, self.inner.rstop)
    }

    /// Trim the read, keeping only the part before the match.
    pub fn trimmed(&self, read: &SequenceRecord) -> SequenceRecord {
        read.slice(0, self.inner.rstart)
    }

    /// The base adjacent to the adapter (just before the match start).
    pub fn adjacent_base(&self) -> &str {
        if self.inner.rstart > 0 {
            &self.inner.sequence[self.inner.rstart - 1..self.inner.rstart]
        } else {
            ""
        }
    }

    /// Length of the removed sequence.
    pub fn removed_sequence_length(&self) -> usize {
        self.inner.sequence.len() - self.inner.rstart
    }

    /// Convenience accessors
    pub fn score(&self) -> i32 {
        self.inner.score
    }
    pub fn errors(&self) -> i32 {
        self.inner.errors
    }
    pub fn rstart(&self) -> usize {
        self.inner.rstart
    }
    pub fn rstop(&self) -> usize {
        self.inner.rstop
    }

    pub fn match_sequence(&self) -> &str {
        self.inner.match_sequence()
    }

    pub fn get_info_records(&self, read: &SequenceRecord) -> Vec<Vec<String>> {
        vec![self.inner.get_info_record(read)]
    }
}

/// A match from a LinkedAdapter, combining a front and a back match.
#[derive(Debug, Clone)]
pub struct LinkedMatch {
    pub front_match: Option<RemoveBeforeMatch>,
    pub back_match: Option<RemoveAfterMatch>,
    pub adapter_name: String,
}

impl LinkedMatch {
    pub fn new(
        front_match: Option<RemoveBeforeMatch>,
        back_match: Option<RemoveAfterMatch>,
        adapter_name: &str,
    ) -> Self {
        assert!(
            front_match.is_some() || back_match.is_some(),
            "At least one of front_match or back_match must be Some"
        );
        Self {
            front_match,
            back_match,
            adapter_name: adapter_name.to_string(),
        }
    }

    pub fn score(&self) -> i32 {
        let mut s = 0;
        if let Some(ref m) = self.front_match {
            s += m.score();
        }
        if let Some(ref m) = self.back_match {
            s += m.score();
        }
        s
    }

    pub fn errors(&self) -> i32 {
        let mut e = 0;
        if let Some(ref m) = self.front_match {
            e += m.errors();
        }
        if let Some(ref m) = self.back_match {
            e += m.errors();
        }
        e
    }

    pub fn trimmed(&self, read: &SequenceRecord) -> SequenceRecord {
        let mut result = read.clone();
        if let Some(ref fm) = self.front_match {
            result = fm.trimmed(&result);
        }
        if let Some(ref bm) = self.back_match {
            result = bm.trimmed(&result);
        }
        result
    }

    pub fn remainder_interval(&self) -> (usize, usize) {
        let mut matches: Vec<(usize, usize)> = Vec::new();
        if let Some(ref m) = self.front_match {
            matches.push(m.remainder_interval());
        }
        if let Some(ref m) = self.back_match {
            matches.push(m.remainder_interval());
        }
        remainder(&matches)
    }

    pub fn retained_adapter_interval(&self) -> (usize, usize) {
        let (start, offset) = if let Some(ref fm) = self.front_match {
            (fm.rstart(), fm.rstop())
        } else {
            (0, 0)
        };
        let end = if let Some(ref bm) = self.back_match {
            bm.rstop() + offset
        } else if let Some(ref fm) = self.front_match {
            fm.inner.sequence.len()
        } else {
            0
        };
        (start, end)
    }

    pub fn get_info_records(&self, read: &SequenceRecord) -> Vec<Vec<String>> {
        let mut records = Vec::new();
        let mut current_read = read.clone();

        if let Some(ref fm) = self.front_match {
            let mut record = fm.inner.get_info_record(read);
            if record.len() > 7 {
                record[7] = format!("{};1", self.adapter_name);
            }
            records.push(record);
            current_read = fm.trimmed(&current_read);
        }
        if let Some(ref bm) = self.back_match {
            let mut record = bm.inner.get_info_record(&current_read);
            if record.len() > 7 {
                record[7] = format!("{};2", self.adapter_name);
            }
            records.push(record);
        }
        records
    }

    pub fn match_sequence(&self) -> String {
        let front_seq = self
            .front_match
            .as_ref()
            .map(|m| m.match_sequence())
            .unwrap_or("");
        let back_seq = self
            .back_match
            .as_ref()
            .map(|m| m.match_sequence())
            .unwrap_or("");
        format!("{},{}", front_seq, back_seq)
    }
}

/// Enum dispatch for all match types.
#[derive(Debug, Clone)]
pub enum AdapterMatch {
    RemoveBefore(RemoveBeforeMatch),
    RemoveAfter(RemoveAfterMatch),
    Linked(LinkedMatch),
}

impl AdapterMatch {
    pub fn score(&self) -> i32 {
        match self {
            AdapterMatch::RemoveBefore(m) => m.score(),
            AdapterMatch::RemoveAfter(m) => m.score(),
            AdapterMatch::Linked(m) => m.score(),
        }
    }

    pub fn errors(&self) -> i32 {
        match self {
            AdapterMatch::RemoveBefore(m) => m.errors(),
            AdapterMatch::RemoveAfter(m) => m.errors(),
            AdapterMatch::Linked(m) => m.errors(),
        }
    }

    pub fn trimmed(&self, read: &SequenceRecord) -> SequenceRecord {
        match self {
            AdapterMatch::RemoveBefore(m) => m.trimmed(read),
            AdapterMatch::RemoveAfter(m) => m.trimmed(read),
            AdapterMatch::Linked(m) => m.trimmed(read),
        }
    }

    pub fn remainder_interval(&self) -> (usize, usize) {
        match self {
            AdapterMatch::RemoveBefore(m) => m.remainder_interval(),
            AdapterMatch::RemoveAfter(m) => m.remainder_interval(),
            AdapterMatch::Linked(m) => m.remainder_interval(),
        }
    }

    pub fn retained_adapter_interval(&self) -> (usize, usize) {
        match self {
            AdapterMatch::RemoveBefore(m) => m.retained_adapter_interval(),
            AdapterMatch::RemoveAfter(m) => m.retained_adapter_interval(),
            AdapterMatch::Linked(m) => m.retained_adapter_interval(),
        }
    }

    pub fn get_info_records(&self, read: &SequenceRecord) -> Vec<Vec<String>> {
        match self {
            AdapterMatch::RemoveBefore(m) => m.get_info_records(read),
            AdapterMatch::RemoveAfter(m) => m.get_info_records(read),
            AdapterMatch::Linked(m) => m.get_info_records(read),
        }
    }

    pub fn match_sequence(&self) -> String {
        match self {
            AdapterMatch::RemoveBefore(m) => m.match_sequence().to_string(),
            AdapterMatch::RemoveAfter(m) => m.match_sequence().to_string(),
            AdapterMatch::Linked(m) => m.match_sequence(),
        }
    }
}

/// Determine the section of the read that remains after applying multiple matches.
fn remainder(intervals: &[(usize, usize)]) -> (usize, usize) {
    assert!(!intervals.is_empty(), "intervals must not be empty");
    let mut start = 0usize;
    let mut last_stop = 0usize;
    let mut last_start = 0usize;
    for &(match_start, match_stop) in intervals {
        start += match_start;
        last_start = match_start;
        last_stop = match_stop;
    }
    let length = last_stop - last_start;
    (start, start + length)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_remove_before_match() {
        let m = RemoveBeforeMatch::new(0, 5, 0, 5, 5, 0, "ACGTAACGTACGT", "adapter1", "ACGTA");
        assert_eq!(m.remainder_interval(), (5, 13));
        assert_eq!(m.removed_sequence_length(), 5);
        assert_eq!(m.rest(), "");
        assert_eq!(m.match_sequence(), "ACGTA");
    }

    #[test]
    fn test_remove_after_match() {
        let m = RemoveAfterMatch::new(0, 5, 8, 13, 5, 0, "ACGTACGTACGTA", "adapter1", "ACGTA");
        assert_eq!(m.remainder_interval(), (0, 8));
        assert_eq!(m.removed_sequence_length(), 5);
        assert_eq!(m.rest(), "");
        assert_eq!(m.adjacent_base(), "T");
    }

    #[test]
    fn test_remove_after_match_adjacent_base_empty() {
        let m = RemoveAfterMatch::new(0, 5, 0, 5, 5, 0, "ACGTA", "adapter1", "ACGTA");
        assert_eq!(m.adjacent_base(), "");
    }

    #[test]
    fn test_linked_match() {
        let fm = RemoveBeforeMatch::new(0, 3, 0, 3, 3, 0, "ACGTACGTACGT", "ad", "ACG");
        let bm = RemoveAfterMatch::new(0, 3, 6, 9, 3, 0, "TACGTACGT", "ad", "CGT");
        let lm = LinkedMatch::new(Some(fm), Some(bm), "linked1");
        assert_eq!(lm.score(), 6);
        assert_eq!(lm.errors(), 0);
    }

    #[test]
    fn test_adapter_match_enum() {
        let m = RemoveBeforeMatch::new(0, 3, 0, 3, 3, 0, "ACGTACGT", "ad", "ACG");
        let am = AdapterMatch::RemoveBefore(m);
        assert_eq!(am.score(), 3);
        assert_eq!(am.errors(), 0);
    }

    #[test]
    fn test_remainder() {
        assert_eq!(remainder(&[(5, 10)]), (5, 10));
        assert_eq!(remainder(&[(3, 10), (0, 5)]), (3, 8));
    }

    #[test]
    fn test_trimmed_remove_before() {
        let m = RemoveBeforeMatch::new(0, 3, 0, 3, 3, 0, "ACGTACGT", "ad", "ACG");
        let read = SequenceRecord::new("r1", "ACGTACGT", Some("IIIIIIII"));
        let trimmed = m.trimmed(&read);
        assert_eq!(trimmed.sequence, "TACGT");
    }

    #[test]
    fn test_trimmed_remove_after() {
        let m = RemoveAfterMatch::new(0, 3, 5, 8, 3, 0, "ACGTACGT", "ad", "CGT");
        let read = SequenceRecord::new("r1", "ACGTACGT", Some("IIIIIIII"));
        let trimmed = m.trimmed(&read);
        assert_eq!(trimmed.sequence, "ACGTA");
    }
}
