// Linked adapter: a 5' adapter combined with a 3' adapter.
//
// Port of LinkedAdapter from cutadapt/adapters.py.

use super::matches::{AdapterMatch, LinkedMatch};
use super::single::SingleAdapter;
use super::statistics::{AdapterStatistics, EndStatistics, LinkedAdapterStatistics};

/// A linked adapter: a front adapter paired with a back adapter.
///
/// The front adapter is matched first. If it matches, the portion of the
/// read after it is searched for the back adapter.
#[derive(Debug)]
pub struct LinkedAdapter {
    pub name: String,
    pub front_adapter: SingleAdapter,
    pub back_adapter: SingleAdapter,
    pub front_required: bool,
    pub back_required: bool,
}

impl LinkedAdapter {
    pub fn new(
        front_adapter: SingleAdapter,
        back_adapter: SingleAdapter,
        front_required: bool,
        back_required: bool,
        name: Option<String>,
    ) -> Self {
        let name = name.unwrap_or_else(|| super::generate_adapter_name());
        Self {
            name,
            front_adapter,
            back_adapter,
            front_required,
            back_required,
        }
    }

    pub fn descriptive_identifier(&self) -> &str {
        "linked"
    }

    pub fn description(&self) -> &str {
        "linked"
    }

    /// The combined sequence representation.
    pub fn sequence(&self) -> String {
        format!(
            "{}...{}",
            self.front_adapter.sequence(),
            self.back_adapter.sequence()
        )
    }

    /// Spec representation.
    pub fn spec(&self) -> String {
        format!(
            "{}...{}",
            self.front_adapter.spec(),
            self.back_adapter.spec()
        )
    }

    /// Attempt to match the linked adapters against a sequence.
    pub fn match_to(&mut self, sequence: &str) -> Option<AdapterMatch> {
        let front_match = self.front_adapter.match_to(sequence);
        if self.front_required && front_match.is_none() {
            return None;
        }

        // If the front adapter matched, search only the remaining part for the back adapter.
        let back_sequence: String;
        let (front_rbm, trimmed_seq) = if let Some(ref fm) = front_match {
            match fm {
                AdapterMatch::RemoveBefore(ref rbm) => {
                    let (start, _) = rbm.remainder_interval();
                    back_sequence = sequence[start..].to_string();
                    (Some(rbm.clone()), back_sequence.as_str())
                }
                _ => {
                    back_sequence = sequence.to_string();
                    (None, back_sequence.as_str())
                }
            }
        } else {
            back_sequence = sequence.to_string();
            (None, back_sequence.as_str())
        };

        let back_match = self.back_adapter.match_to(trimmed_seq);
        let back_ram = back_match.and_then(|bm| match bm {
            AdapterMatch::RemoveAfter(ram) => Some(ram),
            _ => None,
        });

        if back_ram.is_none() && (self.back_required || front_rbm.is_none()) {
            return None;
        }

        Some(AdapterMatch::Linked(LinkedMatch::new(
            front_rbm,
            back_ram,
            &self.name,
        )))
    }

    /// Create statistics for this linked adapter.
    pub fn create_statistics(&self) -> Box<dyn AdapterStatistics> {
        let front_es = EndStatistics::new(
            self.front_adapter.max_error_rate(),
            self.front_adapter.sequence(),
            self.front_adapter.effective_length(),
            self.front_adapter.adapter_wildcards(),
            self.front_adapter.indels(),
            self.front_adapter.descriptive_identifier(),
            self.front_adapter.allows_partial_matches(),
            true, // remove_prefix
        );
        let back_es = EndStatistics::new(
            self.back_adapter.max_error_rate(),
            self.back_adapter.sequence(),
            self.back_adapter.effective_length(),
            self.back_adapter.adapter_wildcards(),
            self.back_adapter.indels(),
            self.back_adapter.descriptive_identifier(),
            self.back_adapter.allows_partial_matches(),
            false, // remove_prefix
        );
        Box::new(LinkedAdapterStatistics::new(&self.name, front_es, back_es))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::single::*;

    fn make_front(seq: &str) -> SingleAdapter {
        let params = SingleAdapterParams {
            sequence: seq.to_string(),
            max_errors: 0.1,
            min_overlap: 3,
            read_wildcards: false,
            adapter_wildcards: true,
            name: Some("front".to_string()),
            indels: true,
        };
        SingleAdapter::Front(FrontAdapter::new(params))
    }

    fn make_back(seq: &str) -> SingleAdapter {
        let params = SingleAdapterParams {
            sequence: seq.to_string(),
            max_errors: 0.1,
            min_overlap: 3,
            read_wildcards: false,
            adapter_wildcards: true,
            name: Some("back".to_string()),
            indels: true,
        };
        SingleAdapter::Back(BackAdapter::new(params))
    }

    #[test]
    fn test_linked_adapter_both_match() {
        let front = make_front("AAAA");
        let back = make_back("TTTT");
        let mut linked = LinkedAdapter::new(front, back, true, true, Some("linked".to_string()));
        // Construct a sequence with front at start and back at end
        let m = linked.match_to("AAAACCCCCTTTT");
        assert!(m.is_some());
        let m = m.unwrap();
        assert!(m.score() > 0);
    }

    #[test]
    fn test_linked_adapter_front_required_missing() {
        let front = make_front("AAAA");
        let back = make_back("TTTT");
        let mut linked = LinkedAdapter::new(front, back, true, true, Some("linked".to_string()));
        let m = linked.match_to("CCCCCTTTT");
        // Front is required but not found -> None
        assert!(m.is_none());
    }

    #[test]
    fn test_linked_adapter_front_optional() {
        let front = make_front("AAAA");
        let back = make_back("TTTT");
        let mut linked = LinkedAdapter::new(front, back, false, true, Some("linked".to_string()));
        let m = linked.match_to("CCCCCTTTT");
        // Front is optional, back found -> match
        assert!(m.is_some());
    }

    #[test]
    fn test_linked_adapter_descriptive() {
        let front = make_front("AAAA");
        let back = make_back("TTTT");
        let linked = LinkedAdapter::new(front, back, true, true, Some("linked".to_string()));
        assert_eq!(linked.descriptive_identifier(), "linked");
        assert_eq!(linked.sequence(), "AAAA...TTTT");
    }

    #[test]
    fn test_linked_adapter_statistics() {
        let front = make_front("AAAA");
        let back = make_back("TTTT");
        let linked = LinkedAdapter::new(front, back, true, true, Some("linked".to_string()));
        let stats = linked.create_statistics();
        assert_eq!(stats.name(), "linked");
        let (front_end, back_end) = stats.end_statistics();
        assert!(front_end.is_some());
        assert!(back_end.is_some());
    }
}
