/// Read-name modification modifiers.
///
/// Partial port of the renaming modifiers from cutadapt/modifiers.py.
/// The full `Renamer` with its template compilation is stubbed for now.

use regex::Regex;

use crate::info::ModificationInfo;
use crate::record::SequenceRecord;

use super::SingleEndModifier;

// ---------------------------------------------------------------------------
// LengthTagModifier
// ---------------------------------------------------------------------------

/// Replace `length=<N>` tags in the read name with the current sequence
/// length.
///
/// Equivalent to Python's `LengthTagModifier`.
pub struct LengthTagModifier {
    pub length_tag: String,
    regex: Regex,
}

impl LengthTagModifier {
    pub fn new(length_tag: &str) -> Self {
        let pattern = format!(r"\b{}{}\b", regex::escape(length_tag), r"[0-9]*");
        Self {
            length_tag: length_tag.to_string(),
            regex: Regex::new(&pattern).expect("Invalid length_tag pattern"),
        }
    }
}

impl SingleEndModifier for LengthTagModifier {
    fn modify(
        &mut self,
        mut read: SequenceRecord,
        _info: &mut ModificationInfo,
    ) -> SequenceRecord {
        if read.name.contains(&self.length_tag) {
            let replacement = format!("{}{}", self.length_tag, read.sequence.len());
            read.name = self.regex.replace_all(&read.name, replacement.as_str()).to_string();
        }
        read
    }
}

// ---------------------------------------------------------------------------
// SuffixRemover
// ---------------------------------------------------------------------------

/// Remove a fixed suffix from the read name, if present.
///
/// Equivalent to Python's `SuffixRemover`.
pub struct SuffixRemover {
    pub suffix: String,
}

impl SuffixRemover {
    pub fn new(suffix: &str) -> Self {
        Self {
            suffix: suffix.to_string(),
        }
    }
}

impl SingleEndModifier for SuffixRemover {
    fn modify(
        &mut self,
        mut read: SequenceRecord,
        _info: &mut ModificationInfo,
    ) -> SequenceRecord {
        if read.name.ends_with(&self.suffix) {
            let new_len = read.name.len() - self.suffix.len();
            read.name.truncate(new_len);
        }
        read
    }
}

// ---------------------------------------------------------------------------
// PrefixSuffixAdder
// ---------------------------------------------------------------------------

/// Add a prefix and/or suffix to the read name.
///
/// The literal `{name}` in prefix/suffix is replaced with the last matched
/// adapter name (or `"no_adapter"` when there are no matches).
///
/// Equivalent to Python's `PrefixSuffixAdder`.
pub struct PrefixSuffixAdder {
    pub prefix: String,
    pub suffix: String,
}

impl PrefixSuffixAdder {
    pub fn new(prefix: &str, suffix: &str) -> Self {
        Self {
            prefix: prefix.to_string(),
            suffix: suffix.to_string(),
        }
    }
}

impl SingleEndModifier for PrefixSuffixAdder {
    fn modify(
        &mut self,
        mut read: SequenceRecord,
        info: &mut ModificationInfo,
    ) -> SequenceRecord {
        // In the full implementation, adapter_name comes from the last match.
        // For now, we use "no_adapter" when there are no matches.
        let adapter_name = if info.matches.is_empty() {
            "no_adapter"
        } else {
            // Placeholder — the real adapter name will come from the Match type
            "adapter"
        };
        let prefix = self.prefix.replace("{name}", adapter_name);
        let suffix = self.suffix.replace("{name}", adapter_name);
        read.name = format!("{}{}{}", prefix, read.name, suffix);
        read
    }
}

// ---------------------------------------------------------------------------
// Renamer (stub)
// ---------------------------------------------------------------------------

/// Rename reads using a template string.
///
/// This is a simplified stub. The full Python implementation compiles a
/// rename function from a template with placeholders like `{header}`,
/// `{id}`, `{comment}`, `{adapter_name}`, etc.
///
/// Currently only performs a plain string replacement of `{header}` with the
/// full read name, `{id}` with the part before the first whitespace, and
/// `{comment}` with the rest.
pub struct Renamer {
    pub template: String,
}

impl Renamer {
    pub fn new(template: &str) -> Self {
        Self {
            template: template.replace(r"\t", "\t"),
        }
    }

    /// Split a read name into (id, comment).
    fn parse_name(name: &str) -> (&str, &str) {
        match name.split_once(char::is_whitespace) {
            Some((id, comment)) => (id, comment),
            None => (name, ""),
        }
    }
}

impl SingleEndModifier for Renamer {
    fn modify(
        &mut self,
        mut read: SequenceRecord,
        info: &mut ModificationInfo,
    ) -> SequenceRecord {
        let (id, comment) = Self::parse_name(&read.name);
        let id = id.to_string();
        let comment = comment.to_string();

        let cut_prefix = info
            .cut_prefix
            .as_deref()
            .unwrap_or("")
            .to_string();
        let cut_suffix = info
            .cut_suffix
            .as_deref()
            .unwrap_or("")
            .to_string();
        let rc = match info.is_rc {
            Some(true) => "rc",
            _ => "",
        };
        let adapter_name = if info.matches.is_empty() {
            "no_adapter"
        } else {
            "adapter"
        };

        let new_name = self
            .template
            .replace("{header}", &read.name)
            .replace("{id}", &id)
            .replace("{comment}", &comment)
            .replace("{cut_prefix}", &cut_prefix)
            .replace("{cut_suffix}", &cut_suffix)
            .replace("{adapter_name}", adapter_name)
            .replace("{rc}", rc);

        read.name = new_name;
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

    // -- LengthTagModifier --

    #[test]
    fn test_length_tag_modifier_replaces() {
        let read = SequenceRecord::new("read1 length=100", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut modifier = LengthTagModifier::new("length=");

        let result = modifier.modify(read, &mut info);
        assert_eq!(result.name, "read1 length=4");
    }

    #[test]
    fn test_length_tag_modifier_no_tag() {
        let read = SequenceRecord::new("read1", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut modifier = LengthTagModifier::new("length=");

        let result = modifier.modify(read, &mut info);
        assert_eq!(result.name, "read1");
    }

    #[test]
    fn test_length_tag_modifier_zero_value() {
        // "length=0" has a digit so the \b after [0-9]* will match
        let read = SequenceRecord::new("read1 length=0", "ACGTAC", Some("IIIIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut modifier = LengthTagModifier::new("length=");

        let result = modifier.modify(read, &mut info);
        assert_eq!(result.name, "read1 length=6");
    }

    // -- SuffixRemover --

    #[test]
    fn test_suffix_remover_removes() {
        let read = SequenceRecord::new("read1/1", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut modifier = SuffixRemover::new("/1");

        let result = modifier.modify(read, &mut info);
        assert_eq!(result.name, "read1");
    }

    #[test]
    fn test_suffix_remover_no_match() {
        let read = SequenceRecord::new("read1", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut modifier = SuffixRemover::new("/1");

        let result = modifier.modify(read, &mut info);
        assert_eq!(result.name, "read1");
    }

    // -- PrefixSuffixAdder --

    #[test]
    fn test_prefix_suffix_adder() {
        let read = SequenceRecord::new("read1", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut modifier = PrefixSuffixAdder::new("pre_", "_suf");

        let result = modifier.modify(read, &mut info);
        assert_eq!(result.name, "pre_read1_suf");
    }

    #[test]
    fn test_prefix_suffix_adder_with_name_placeholder() {
        let read = SequenceRecord::new("read1", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        // No matches, so {name} becomes "no_adapter"
        let mut modifier = PrefixSuffixAdder::new("{name}_", "_{name}");

        let result = modifier.modify(read, &mut info);
        assert_eq!(result.name, "no_adapter_read1_no_adapter");
    }

    // -- Renamer --

    #[test]
    fn test_renamer_header() {
        let read = SequenceRecord::new("read1 some comment", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut renamer = Renamer::new("{header}");

        let result = renamer.modify(read, &mut info);
        assert_eq!(result.name, "read1 some comment");
    }

    #[test]
    fn test_renamer_id_comment() {
        let read = SequenceRecord::new("read1 some comment", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut renamer = Renamer::new("{id}_{comment}");

        let result = renamer.modify(read, &mut info);
        assert_eq!(result.name, "read1_some comment");
    }

    #[test]
    fn test_renamer_id_only() {
        let read = SequenceRecord::new("read1", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut renamer = Renamer::new("{id}");

        let result = renamer.modify(read, &mut info);
        assert_eq!(result.name, "read1");
    }

    #[test]
    fn test_renamer_cut_prefix_suffix() {
        let read = SequenceRecord::new("read1", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        info.cut_prefix = Some("AC".to_string());
        info.cut_suffix = Some("GT".to_string());
        let mut renamer = Renamer::new("{id}_p{cut_prefix}_s{cut_suffix}");

        let result = renamer.modify(read, &mut info);
        assert_eq!(result.name, "read1_pAC_sGT");
    }

    #[test]
    fn test_renamer_rc() {
        let read = SequenceRecord::new("read1", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        info.is_rc = Some(true);
        let mut renamer = Renamer::new("{id}_{rc}");

        let result = renamer.modify(read, &mut info);
        assert_eq!(result.name, "read1_rc");
    }

    #[test]
    fn test_renamer_rc_false() {
        let read = SequenceRecord::new("read1", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        info.is_rc = Some(false);
        let mut renamer = Renamer::new("{id}_{rc}");

        let result = renamer.modify(read, &mut info);
        assert_eq!(result.name, "read1_");
    }

    #[test]
    fn test_renamer_tab_escape() {
        let read = SequenceRecord::new("read1", "ACGT", Some("IIII"));
        let mut info = ModificationInfo::new(read.clone());
        let mut renamer = Renamer::new("{id}\\t{comment}");
        // The \\t in the raw string should have been converted to \t during new()
        // But since the input is a Rust string literal, the backslash-t is literal chars \ and t.
        // The Renamer::new replaces r"\t" with "\t".
        let renamer2 = Renamer::new(r"{id}\t{comment}");
        assert!(renamer2.template.contains('\t'));

        let result = renamer.modify(read, &mut info);
        // In this case the template was "{id}\\t{comment}" which after Rust string
        // processing is "{id}\t{comment}" which Renamer::new converts to
        // "{id}<TAB>{comment}".
        assert!(result.name.contains('\t'));
    }
}
