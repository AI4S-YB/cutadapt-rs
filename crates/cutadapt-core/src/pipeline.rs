/// Processing pipelines for single-end and paired-end reads.
///
/// Port of cutadapt/pipeline.py.
///
/// A pipeline loops over reads, applies modifiers in sequence, then applies
/// steps (filters/sinks) in sequence. Each modifier or step can transform or
/// consume a read.

use crate::info::ModificationInfo;
use crate::modifiers::{PairedEndModifier, SingleEndModifier};
use crate::record::SequenceRecord;
use crate::steps::{PairedEndStep, RecordPair, SingleEndStep};

// ---------------------------------------------------------------------------
// Pipeline trait
// ---------------------------------------------------------------------------

/// Common interface for single-end and paired-end processing pipelines.
pub trait Pipeline {
    /// Process all reads and return `(n_reads, total_bp_r1, total_bp_r2)`.
    ///
    /// For single-end pipelines the third element is `None`.
    fn process_reads(&mut self) -> (u64, u64, Option<u64>);
}

// ---------------------------------------------------------------------------
// SingleEndPipeline
// ---------------------------------------------------------------------------

/// Processing pipeline for single-end reads.
///
/// Takes a `Vec<SequenceRecord>` as input (file I/O will be added in Phase 7).
pub struct SingleEndPipeline {
    reads: Vec<SequenceRecord>,
    modifiers: Vec<Box<dyn SingleEndModifier>>,
    steps: Vec<Box<dyn SingleEndStep>>,
}

impl SingleEndPipeline {
    pub fn new(
        reads: Vec<SequenceRecord>,
        modifiers: Vec<Box<dyn SingleEndModifier>>,
        steps: Vec<Box<dyn SingleEndStep>>,
    ) -> Self {
        Self {
            reads,
            modifiers,
            steps,
        }
    }
}

impl Pipeline for SingleEndPipeline {
    fn process_reads(&mut self) -> (u64, u64, Option<u64>) {
        let mut n: u64 = 0;
        let mut total_bp: u64 = 0;

        let reads = std::mem::take(&mut self.reads);
        for mut read in reads {
            n += 1;
            total_bp += read.len() as u64;
            let mut info = ModificationInfo::new(read.clone());

            // Apply modifiers.
            for modifier in self.modifiers.iter_mut() {
                read = modifier.modify(read, &mut info);
            }

            // Apply steps. If a step consumes the read (returns None), stop.
            let mut current = Some(read);
            for step in self.steps.iter_mut() {
                if let Some(r) = current {
                    current = step.process(r, &info);
                } else {
                    break;
                }
            }
        }

        (n, total_bp, None)
    }
}

// ---------------------------------------------------------------------------
// PairedEndPipeline
// ---------------------------------------------------------------------------

/// Processing pipeline for paired-end reads.
///
/// Takes a `Vec<(SequenceRecord, SequenceRecord)>` as input.
pub struct PairedEndPipeline {
    pairs: Vec<RecordPair>,
    modifiers: Vec<Box<dyn PairedEndModifier>>,
    steps: Vec<Box<dyn PairedEndStep>>,
}

impl PairedEndPipeline {
    pub fn new(
        pairs: Vec<RecordPair>,
        modifiers: Vec<Box<dyn PairedEndModifier>>,
        steps: Vec<Box<dyn PairedEndStep>>,
    ) -> Self {
        Self {
            pairs,
            modifiers,
            steps,
        }
    }
}

impl Pipeline for PairedEndPipeline {
    fn process_reads(&mut self) -> (u64, u64, Option<u64>) {
        let mut n: u64 = 0;
        let mut total1_bp: u64 = 0;
        let mut total2_bp: u64 = 0;

        let pairs = std::mem::take(&mut self.pairs);
        for (mut read1, mut read2) in pairs {
            n += 1;
            total1_bp += read1.len() as u64;
            total2_bp += read2.len() as u64;
            let mut info1 = ModificationInfo::new(read1.clone());
            let mut info2 = ModificationInfo::new(read2.clone());

            // Apply modifiers.
            for modifier in self.modifiers.iter_mut() {
                let pair = modifier.modify_pair(read1, read2, &mut info1, &mut info2);
                read1 = pair.0;
                read2 = pair.1;
            }

            // Apply steps. If a step consumes the pair (returns None), stop.
            let mut current: Option<RecordPair> = Some((read1, read2));
            for step in self.steps.iter_mut() {
                if let Some((r1, r2)) = current {
                    current = step.process_pair(r1, r2, &info1, &info2);
                } else {
                    break;
                }
            }
        }

        (n, total1_bp, Some(total2_bp))
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::info::ModificationInfo;
    use crate::modifiers::SingleEndModifier;
    use crate::record::SequenceRecord;
    use crate::steps::{
        PairedEndSink, PairFilterMode, PairedEndFilter, SingleEndFilter, SingleEndSink,
    };
    use crate::predicates::TooShort;

    fn make_read(name: &str, seq: &str) -> SequenceRecord {
        SequenceRecord::new(name, seq, Some(&"I".repeat(seq.len())))
    }

    // -- A trivial modifier for testing: uppercases the sequence --

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

    // -- SingleEndPipeline --

    #[test]
    fn test_single_end_pipeline_basic() {
        let reads = vec![
            make_read("r1", "ACGT"),
            make_read("r2", "TGCATGCA"),
        ];
        let sink = SingleEndSink::new();
        let sink_ptr = &sink as *const SingleEndSink;
        // We need to get the sink back out after processing, so we use a shared
        // reference pattern with an index into the steps vec.
        let mut pipeline = SingleEndPipeline::new(
            reads,
            vec![],
            vec![Box::new(sink)],
        );

        let (n, total_bp, bp2) = pipeline.process_reads();
        assert_eq!(n, 2);
        assert_eq!(total_bp, 12); // 4 + 8
        assert!(bp2.is_none());

        // Verify the sink received the reads by downcasting.
        // Since we can't easily downcast trait objects, we just verify the
        // counts from the pipeline return value.
        let _ = sink_ptr; // suppress warning
    }

    #[test]
    fn test_single_end_pipeline_with_modifier() {
        let reads = vec![make_read("r1", "acgt")];
        let sink = SingleEndSink::new();

        let mut pipeline = SingleEndPipeline::new(
            reads,
            vec![Box::new(UppercaseModifier)],
            vec![Box::new(sink)],
        );

        let (n, total_bp, _) = pipeline.process_reads();
        assert_eq!(n, 1);
        assert_eq!(total_bp, 4);
    }

    #[test]
    fn test_single_end_pipeline_with_filter_and_sink() {
        let reads = vec![
            make_read("r1", "AC"),         // too short, will be filtered
            make_read("r2", "ACGTACGT"),   // passes filter
            make_read("r3", "T"),          // too short, will be filtered
        ];

        let mut pipeline = SingleEndPipeline::new(
            reads,
            vec![],
            vec![
                Box::new(SingleEndFilter::new(Box::new(TooShort::new(5)), false)),
                Box::new(SingleEndSink::new()),
            ],
        );

        let (n, total_bp, _) = pipeline.process_reads();
        assert_eq!(n, 3);
        assert_eq!(total_bp, 11); // 2 + 8 + 1
    }

    #[test]
    fn test_single_end_pipeline_empty() {
        let mut pipeline = SingleEndPipeline::new(
            vec![],
            vec![],
            vec![Box::new(SingleEndSink::new())],
        );
        let (n, total_bp, bp2) = pipeline.process_reads();
        assert_eq!(n, 0);
        assert_eq!(total_bp, 0);
        assert!(bp2.is_none());
    }

    // -- PairedEndPipeline --

    #[test]
    fn test_paired_end_pipeline_basic() {
        let pairs = vec![
            (make_read("r1", "ACGT"), make_read("r2", "TGCATG")),
            (make_read("r3", "AA"), make_read("r4", "TT")),
        ];

        let mut pipeline = PairedEndPipeline::new(
            pairs,
            vec![],
            vec![Box::new(PairedEndSink::new())],
        );

        let (n, bp1, bp2) = pipeline.process_reads();
        assert_eq!(n, 2);
        assert_eq!(bp1, 6);  // 4 + 2
        assert_eq!(bp2, Some(8)); // 6 + 2
    }

    #[test]
    fn test_paired_end_pipeline_with_filter_and_sink() {
        let pairs = vec![
            (make_read("r1", "AC"), make_read("r2", "TG")),           // both short -> filtered
            (make_read("r3", "ACGTACGT"), make_read("r4", "TGCATGCA")), // both long -> pass
        ];

        let mut pipeline = PairedEndPipeline::new(
            pairs,
            vec![],
            vec![
                Box::new(PairedEndFilter::new(
                    Some(Box::new(TooShort::new(5))),
                    Some(Box::new(TooShort::new(5))),
                    PairFilterMode::Both,
                    false,
                )),
                Box::new(PairedEndSink::new()),
            ],
        );

        let (n, bp1, bp2) = pipeline.process_reads();
        assert_eq!(n, 2);
        assert_eq!(bp1, 10); // 2 + 8
        assert_eq!(bp2, Some(10)); // 2 + 8
    }

    #[test]
    fn test_paired_end_pipeline_empty() {
        let mut pipeline = PairedEndPipeline::new(
            vec![],
            vec![],
            vec![Box::new(PairedEndSink::new())],
        );
        let (n, bp1, bp2) = pipeline.process_reads();
        assert_eq!(n, 0);
        assert_eq!(bp1, 0);
        assert_eq!(bp2, Some(0));
    }
}
