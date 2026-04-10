/// Adapter trimming modifier.
///
/// Wraps a set of adapters and applies them to a read, trimming the best match.
/// Supports multiple passes (`times > 1`).

use crate::adapters::{Adapter, MultipleAdapters};
use crate::info::{MatchInfo, ModificationInfo};
use crate::record::SequenceRecord;

use super::SingleEndModifier;

/// Applies adapter trimming to a single read.
///
/// On each call to `modify`, up to `times` adapter matches are found and
/// trimmed. Statistics are accumulated across all calls.
pub struct AdapterTrimmer {
    pub adapters: MultipleAdapters,
    pub times: usize,
    /// Total number of reads that had at least one adapter trimmed.
    pub reads_with_adapter: u64,
}

impl AdapterTrimmer {
    pub fn new(adapters: Vec<Adapter>, times: usize) -> Self {
        Self {
            adapters: MultipleAdapters::new(adapters),
            times,
            reads_with_adapter: 0,
        }
    }
}

impl SingleEndModifier for AdapterTrimmer {
    fn modify(
        &mut self,
        mut read: SequenceRecord,
        info: &mut ModificationInfo,
    ) -> SequenceRecord {
        let mut trimmed_any = false;
        for _ in 0..self.times {
            match self.adapters.match_to(&read.sequence) {
                Some((m, idx)) => {
                    let adapter_name = self.adapters.adapters[idx].name().to_string();
                    read = m.trimmed(&read);
                    info.matches.push(MatchInfo {
                        adapter_name: Some(adapter_name),
                    });
                    trimmed_any = true;
                    // If the read is now empty, stop trimming.
                    if read.is_empty() {
                        break;
                    }
                }
                None => break,
            }
        }
        if trimmed_any {
            self.reads_with_adapter += 1;
        }
        read
    }
}
