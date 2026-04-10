/// Per-read pipeline state container.
///
/// Port of cutadapt/info.pyx.

use crate::record::SequenceRecord;

/// An object of this class is created for each read that passes through the pipeline.
/// Any information (except the read itself) that needs to be passed from one modifier
/// to one later in the pipeline or from one modifier to the filters is recorded here.
#[derive(Debug, Clone)]
pub struct ModificationInfo {
    /// List of adapter matches found on this read.
    pub matches: Vec<MatchInfo>,
    /// The original unmodified read.
    pub original_read: SequenceRecord,
    /// Length of prefix cut by UnconditionalCutter.
    pub cut_prefix: Option<String>,
    /// Length of suffix cut by UnconditionalCutter.
    pub cut_suffix: Option<String>,
    /// Whether this read was reverse complemented by ReverseComplementer.
    pub is_rc: Option<bool>,
}

/// Placeholder for match information stored during pipeline processing.
/// Will be replaced with the actual Match type from the adapters module.
#[derive(Debug, Clone)]
pub struct MatchInfo {
    /// Name of the adapter that matched (used for demultiplexing).
    pub adapter_name: Option<String>,
}

impl ModificationInfo {
    pub fn new(read: SequenceRecord) -> Self {
        Self {
            matches: Vec::new(),
            original_read: read,
            cut_prefix: None,
            cut_suffix: None,
            is_rc: None,
        }
    }
}
