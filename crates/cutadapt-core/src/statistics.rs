/// Read length statistics tracking.
///
/// Port of cutadapt/statistics.py.

use std::collections::HashMap;

/// Keep track of the lengths of written reads or read pairs.
#[derive(Debug, Clone, Default)]
pub struct ReadLengthStatistics {
    written_lengths1: HashMap<usize, u64>,
    written_lengths2: HashMap<usize, u64>,
}

impl ReadLengthStatistics {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a single-end read to the statistics.
    pub fn update(&mut self, read_length: usize) {
        *self.written_lengths1.entry(read_length).or_insert(0) += 1;
    }

    /// Add a paired-end read to the statistics.
    pub fn update2(&mut self, read1_length: usize, read2_length: usize) {
        *self.written_lengths1.entry(read1_length).or_insert(0) += 1;
        *self.written_lengths2.entry(read2_length).or_insert(0) += 1;
    }

    /// Return number of written reads or read pairs.
    pub fn written_reads(&self) -> u64 {
        self.written_lengths1.values().sum()
    }

    /// Return total base pairs written as (bp1, bp2).
    pub fn written_bp(&self) -> (u64, u64) {
        (
            Self::compute_total_bp(&self.written_lengths1),
            Self::compute_total_bp(&self.written_lengths2),
        )
    }

    /// Return copies of the length distribution maps.
    pub fn written_lengths(&self) -> (&HashMap<usize, u64>, &HashMap<usize, u64>) {
        (&self.written_lengths1, &self.written_lengths2)
    }

    fn compute_total_bp(counts: &HashMap<usize, u64>) -> u64 {
        counts
            .iter()
            .map(|(&length, &count)| length as u64 * count)
            .sum()
    }
}

impl std::ops::AddAssign for ReadLengthStatistics {
    fn add_assign(&mut self, other: Self) {
        for (length, count) in &other.written_lengths1 {
            *self.written_lengths1.entry(*length).or_insert(0) += count;
        }
        for (length, count) in &other.written_lengths2 {
            *self.written_lengths2.entry(*length).or_insert(0) += count;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_end() {
        let mut stats = ReadLengthStatistics::new();
        stats.update(100);
        stats.update(100);
        stats.update(50);
        assert_eq!(stats.written_reads(), 3);
        assert_eq!(stats.written_bp(), (250, 0));
    }

    #[test]
    fn test_paired_end() {
        let mut stats = ReadLengthStatistics::new();
        stats.update2(100, 90);
        stats.update2(100, 80);
        assert_eq!(stats.written_reads(), 2);
        assert_eq!(stats.written_bp(), (200, 170));
    }

    #[test]
    fn test_add_assign() {
        let mut stats1 = ReadLengthStatistics::new();
        stats1.update(100);
        stats1.update(50);

        let mut stats2 = ReadLengthStatistics::new();
        stats2.update(100);
        stats2.update(75);

        stats1 += stats2;
        assert_eq!(stats1.written_reads(), 4);
        assert_eq!(stats1.written_bp(), (325, 0));
    }

    #[test]
    fn test_empty() {
        let stats = ReadLengthStatistics::new();
        assert_eq!(stats.written_reads(), 0);
        assert_eq!(stats.written_bp(), (0, 0));
    }
}
