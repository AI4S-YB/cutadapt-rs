// Parser for adapter specification strings.
//
// Port of cutadapt/parser.py.
//
// This module parses adapter specification strings like "NAME=SEQUENCE;param=value"
// and creates the appropriate adapter objects.

use std::collections::HashMap;
use std::fmt;
use std::io::BufRead;
use std::path::Path;

use crate::adapters::linked::LinkedAdapter;
use crate::adapters::single::{
    AnywhereAdapter, BackAdapter, FrontAdapter, NonInternalBackAdapter, NonInternalFrontAdapter,
    PrefixAdapter, RightmostBackAdapter, RightmostFrontAdapter, SingleAdapter,
    SingleAdapterParams, SuffixAdapter,
};
use crate::adapters::Adapter;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Errors that can occur during adapter specification parsing.
#[derive(Debug, thiserror::Error)]
pub enum ParseError {
    #[error("{0}")]
    InvalidSpec(String),

    #[error("Unknown parameter '{0}'")]
    UnknownParameter(String),

    #[error("No value given for key '{0}'")]
    MissingValue(String),

    #[error("Key '{0}' specified twice")]
    DuplicateKey(String),

    #[error("'optional' and 'required' cannot be specified at the same time")]
    OptionalAndRequired,

    #[error("'indels' and 'noindels' cannot be specified at the same time")]
    IndelsAndNoindels,

    #[error("{0}")]
    BraceExpansion(String),

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
}

// ---------------------------------------------------------------------------
// SearchParameterValue
// ---------------------------------------------------------------------------

/// A value parsed from a search parameter string. Can be a bool, int, or float.
#[derive(Debug, Clone, PartialEq)]
pub enum ParamValue {
    Bool(bool),
    Int(i64),
    Float(f64),
}

impl ParamValue {
    pub fn as_bool(&self) -> Option<bool> {
        match self {
            ParamValue::Bool(b) => Some(*b),
            _ => None,
        }
    }

    pub fn as_int(&self) -> Option<i64> {
        match self {
            ParamValue::Int(i) => Some(*i),
            _ => None,
        }
    }

    pub fn as_float(&self) -> Option<f64> {
        match self {
            ParamValue::Float(f) => Some(*f),
            ParamValue::Int(i) => Some(*i as f64),
            _ => None,
        }
    }
}

// ---------------------------------------------------------------------------
// parse_search_parameters
// ---------------------------------------------------------------------------

/// Parse `key=value;key=value;...` format into a HashMap.
///
/// Allowed keys and their canonical forms:
/// - `e`, `error_rate`, `max_error_rate` -> `max_errors`
/// - `o` -> `min_overlap`
/// - `max_errors`, `min_overlap`, `anywhere`, `required`, `optional`, `indels`, `noindels`, `rightmost`
pub fn parse_search_parameters(spec: &str) -> Result<HashMap<String, ParamValue>, ParseError> {
    let allowed: HashMap<&str, Option<&str>> = [
        // abbreviations
        ("e", Some("max_errors")),
        ("error_rate", Some("max_errors")),
        ("max_error_rate", Some("max_errors")),
        ("o", Some("min_overlap")),
        // canonical keys
        ("max_errors", None),
        ("min_overlap", None),
        ("anywhere", None),
        ("required", None),
        ("optional", None),
        ("indels", None),
        ("noindels", None),
        ("rightmost", None),
    ]
    .iter()
    .cloned()
    .collect();

    let mut result: HashMap<String, ParamValue> = HashMap::new();

    for field in spec.split(';') {
        let field = field.trim();
        if field.is_empty() {
            continue;
        }
        let (key, value_str) = if let Some(eq_pos) = field.find('=') {
            let k = field[..eq_pos].trim();
            let v = field[eq_pos + 1..].trim();
            if v.is_empty() {
                return Err(ParseError::MissingValue(k.to_string()));
            }
            (k, Some(v))
        } else {
            (field, None)
        };

        if !allowed.contains_key(key) {
            return Err(ParseError::UnknownParameter(key.to_string()));
        }

        // Resolve abbreviations
        let mut canonical = key.to_string();
        while let Some(Some(target)) = allowed.get(canonical.as_str()) {
            canonical = target.to_string();
        }

        let value = if let Some(v) = value_str {
            if let Ok(i) = v.parse::<i64>() {
                ParamValue::Int(i)
            } else if let Ok(f) = v.parse::<f64>() {
                ParamValue::Float(f)
            } else {
                return Err(ParseError::InvalidSpec(format!(
                    "Cannot parse value '{}' for key '{}'",
                    v, key
                )));
            }
        } else {
            ParamValue::Bool(true)
        };

        if result.contains_key(&canonical) {
            return Err(ParseError::DuplicateKey(canonical));
        }
        result.insert(canonical, value);
    }

    if result.contains_key("optional") && result.contains_key("required") {
        return Err(ParseError::OptionalAndRequired);
    }
    if result.contains_key("indels") && result.contains_key("noindels") {
        return Err(ParseError::IndelsAndNoindels);
    }

    // Normalize: optional -> required=false
    if result.remove("optional").is_some() {
        result.insert("required".to_string(), ParamValue::Bool(false));
    }
    // Normalize: noindels -> indels=false
    if result.remove("noindels").is_some() {
        result.insert("indels".to_string(), ParamValue::Bool(false));
    }

    Ok(result)
}

// ---------------------------------------------------------------------------
// expand_braces
// ---------------------------------------------------------------------------

/// Replace all occurrences of `x{n}` (where x is any character) with n
/// occurrences of x. Raise an error if the expression cannot be parsed.
///
/// ```
/// use cutadapt_core::parser::expand_braces;
/// assert_eq!(expand_braces("TGA{5}CT").unwrap(), "TGAAAAACT");
/// ```
pub fn expand_braces(sequence: &str) -> Result<String, ParseError> {
    // Split on { and }, keeping delimiters
    let parts = split_braces(sequence);
    let mut result = String::new();

    #[derive(Debug)]
    enum State {
        Normal,
        AfterOpen,
        GotCount(usize),
    }

    let mut state = State::Normal;

    for s in &parts {
        if s.is_empty() {
            continue;
        }
        match state {
            State::Normal => {
                if *s == "{" {
                    if result.is_empty() {
                        return Err(ParseError::BraceExpansion(
                            "\"{\" must be used after a character".to_string(),
                        ));
                    }
                    state = State::AfterOpen;
                } else if *s == "}" {
                    return Err(ParseError::BraceExpansion(
                        "\"}\" cannot be used here".to_string(),
                    ));
                } else {
                    result.push_str(s);
                }
            }
            State::AfterOpen => {
                let count: usize = s.parse().map_err(|_| {
                    ParseError::BraceExpansion(format!("Expected a number inside braces, got '{}'", s))
                })?;
                if count > 10000 {
                    return Err(ParseError::BraceExpansion(format!(
                        "Value {} invalid",
                        count
                    )));
                }
                state = State::GotCount(count);
            }
            State::GotCount(count) => {
                if *s != "}" {
                    return Err(ParseError::BraceExpansion(
                        "\"}\" expected".to_string(),
                    ));
                }
                // Replace the last character with `count` copies of itself
                if let Some(last_char) = result.pop() {
                    for _ in 0..count {
                        result.push(last_char);
                    }
                }
                state = State::Normal;
            }
        }
    }

    match state {
        State::Normal => Ok(result),
        _ => Err(ParseError::BraceExpansion(
            "Unterminated expression".to_string(),
        )),
    }
}

/// Split a string on `{` and `}`, keeping delimiters as separate items.
fn split_braces(s: &str) -> Vec<&str> {
    let mut result = Vec::new();
    let mut start = 0;
    for (i, c) in s.char_indices() {
        if c == '{' || c == '}' {
            if i > start {
                result.push(&s[start..i]);
            }
            result.push(&s[i..i + 1]);
            start = i + 1;
        }
    }
    if start < s.len() {
        result.push(&s[start..]);
    }
    result
}

// ---------------------------------------------------------------------------
// Restriction
// ---------------------------------------------------------------------------

/// Placement restriction for an adapter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Restriction {
    Anchored,
    NonInternal,
}

impl fmt::Display for Restriction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Restriction::Anchored => write!(f, "anchored"),
            Restriction::NonInternal => write!(f, "noninternal"),
        }
    }
}

// ---------------------------------------------------------------------------
// AdapterType
// ---------------------------------------------------------------------------

/// The type of adapter based on the command-line flag used.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdapterType {
    /// 5' adapter (`-g`)
    Front,
    /// 3' adapter (`-a`)
    Back,
    /// Anywhere adapter (`-b`)
    Anywhere,
}

impl fmt::Display for AdapterType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AdapterType::Front => write!(f, "front"),
            AdapterType::Back => write!(f, "back"),
            AdapterType::Anywhere => write!(f, "anywhere"),
        }
    }
}

impl AdapterType {
    pub fn from_str(s: &str) -> Result<Self, ParseError> {
        match s {
            "front" => Ok(AdapterType::Front),
            "back" => Ok(AdapterType::Back),
            "anywhere" => Ok(AdapterType::Anywhere),
            _ => Err(ParseError::InvalidSpec(format!(
                "adapter_type must be front, back or anywhere, got '{}'",
                s
            ))),
        }
    }
}

// ---------------------------------------------------------------------------
// AdapterSpecification
// ---------------------------------------------------------------------------

/// Description of a single non-linked adapter parsed from a specification string.
#[derive(Debug, Clone, PartialEq)]
pub struct AdapterSpecification {
    pub name: Option<String>,
    pub restriction: Option<Restriction>,
    pub sequence: String,
    pub parameters: HashMap<String, ParamValue>,
    pub adapter_type: AdapterType,
    pub rightmost: bool,
}

impl AdapterSpecification {
    /// Parse an adapter specification string for a non-linked adapter (without `...`)
    /// and return an `AdapterSpecification`.
    ///
    /// Accepted forms:
    /// - `ADAPTER` (basic)
    /// - `NAME=ADAPTER` (named)
    /// - `ADAPTER;param=value` (with parameters)
    /// - `^ADAPTER` or `ADAPTER$` (anchored)
    /// - `XADAPTER` or `ADAPTERX` (non-internal)
    pub fn parse(spec: &str, adapter_type: AdapterType) -> Result<Self, ParseError> {
        let (before_params, parameters_spec) = split_first_semicolon(spec);
        let (name, sequence_part) = extract_name(before_params);
        let sequence_part = sequence_part.trim();
        let mut parameters = parse_search_parameters(parameters_spec)?;
        let sequence_expanded = expand_braces(sequence_part)?;
        let rightmost = parameters.remove("rightmost").is_some();

        // Special case: adapter consisting of only X characters
        if !sequence_expanded.is_empty()
            && sequence_expanded.chars().all(|c| c == 'X' || c == 'x')
        {
            return Ok(AdapterSpecification {
                name,
                restriction: None,
                sequence: sequence_expanded,
                parameters: HashMap::new(),
                adapter_type,
                rightmost: false,
            });
        }

        let (front_restriction, back_restriction, stripped_seq) =
            parse_restrictions(&sequence_expanded)?;

        if adapter_type == AdapterType::Front && back_restriction.is_some() {
            return Err(ParseError::InvalidSpec(
                "Allowed placement restrictions for a 5' adapter are XADAPTER and ^ADAPTER"
                    .to_string(),
            ));
        }
        if adapter_type == AdapterType::Back && front_restriction.is_some() {
            return Err(ParseError::InvalidSpec(
                "Allowed placement restrictions for a 3' adapter are ADAPTERX and ADAPTER$"
                    .to_string(),
            ));
        }

        let restriction = front_restriction.or(back_restriction);

        if adapter_type == AdapterType::Anywhere && restriction.is_some() {
            return Err(ParseError::InvalidSpec(
                "Placement restrictions (with X, ^, $) not supported for 'anywhere' (-b) adapters"
                    .to_string(),
            ));
        }

        if parameters.contains_key("min_overlap") && restriction == Some(Restriction::Anchored) {
            return Err(ParseError::InvalidSpec(
                "Setting 'min_overlap=' (or 'o=') for anchored adapters is not possible because \
                 anchored adapters always need to match in full."
                    .to_string(),
            ));
        }

        // Clamp min_overlap to sequence length
        if let Some(ParamValue::Int(overlap)) = parameters.get("min_overlap") {
            if *overlap as usize > stripped_seq.len() {
                parameters.insert(
                    "min_overlap".to_string(),
                    ParamValue::Int(stripped_seq.len() as i64),
                );
            }
        }

        if rightmost
            && (adapter_type == AdapterType::Anywhere || restriction.is_some())
        {
            return Err(ParseError::InvalidSpec(
                "'rightmost' only allowed with regular 5' and 3' adapters".to_string(),
            ));
        }

        Ok(AdapterSpecification {
            name,
            restriction,
            sequence: stripped_seq,
            parameters,
            adapter_type,
            rightmost,
        })
    }

    /// Return the appropriate `SingleAdapter` variant for this specification.
    pub fn adapter_class(&self) -> AdapterClass {
        restriction_to_class(self.adapter_type, self.restriction, self.rightmost)
    }
}

/// Identifies which concrete adapter type should be created.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdapterClass {
    Front,
    Back,
    Anywhere,
    Prefix,
    Suffix,
    NonInternalFront,
    NonInternalBack,
    RightmostFront,
    RightmostBack,
}

fn restriction_to_class(
    adapter_type: AdapterType,
    restriction: Option<Restriction>,
    rightmost: bool,
) -> AdapterClass {
    match adapter_type {
        AdapterType::Front => {
            if rightmost {
                AdapterClass::RightmostFront
            } else {
                match restriction {
                    None => AdapterClass::Front,
                    Some(Restriction::Anchored) => AdapterClass::Prefix,
                    Some(Restriction::NonInternal) => AdapterClass::NonInternalFront,
                }
            }
        }
        AdapterType::Back => {
            if rightmost {
                AdapterClass::RightmostBack
            } else {
                match restriction {
                    None => AdapterClass::Back,
                    Some(Restriction::Anchored) => AdapterClass::Suffix,
                    Some(Restriction::NonInternal) => AdapterClass::NonInternalBack,
                }
            }
        }
        AdapterType::Anywhere => AdapterClass::Anywhere,
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Split on the first semicolon. Returns (before, after).
/// If no semicolon, `after` is an empty string.
fn split_first_semicolon(s: &str) -> (&str, &str) {
    if let Some(pos) = s.find(';') {
        (&s[..pos], &s[pos + 1..])
    } else {
        (s, "")
    }
}

/// Extract the adapter name from `NAME=SEQUENCE`, returning `(name, sequence)`.
fn extract_name(spec: &str) -> (Option<String>, &str) {
    if let Some(eq_pos) = spec.find('=') {
        let name = spec[..eq_pos].trim();
        let rest = spec[eq_pos + 1..].trim();
        (Some(name.to_string()), rest)
    } else {
        (None, spec)
    }
}

/// Parse the restriction notation from a sequence.
///
/// Returns `(front_restriction, back_restriction, stripped_sequence)`.
fn parse_restrictions(spec: &str) -> Result<(Option<Restriction>, Option<Restriction>, String), ParseError> {
    let mut s = spec.to_string();
    let mut front_restriction: Option<Restriction> = None;

    // Front restrictions: ^ or leading X
    if s.starts_with('^') {
        front_restriction = Some(Restriction::Anchored);
        s = s[1..].to_string();
    }
    if s.starts_with('X') || s.starts_with('x') {
        if front_restriction.is_some() {
            return Err(ParseError::InvalidSpec(
                "You cannot use multiple placement restrictions for an adapter at the same time. \
                 Choose one of ^ADAPTER, ADAPTER$, XADAPTER or ADAPTERX"
                    .to_string(),
            ));
        }
        front_restriction = Some(Restriction::NonInternal);
        s = s.trim_start_matches(|c: char| c == 'X' || c == 'x').to_string();
    }

    let mut back_restriction: Option<Restriction> = None;

    // Back restrictions: $ or trailing X
    if s.ends_with('$') {
        back_restriction = Some(Restriction::Anchored);
        s = s[..s.len() - 1].to_string();
    }
    if s.ends_with('X') || s.ends_with('x') {
        if back_restriction.is_some() {
            return Err(ParseError::InvalidSpec(
                "You cannot use multiple placement restrictions for an adapter at the same time. \
                 Choose one of ^ADAPTER, ADAPTER$, XADAPTER or ADAPTERX"
                    .to_string(),
            ));
        }
        back_restriction = Some(Restriction::NonInternal);
        s = s.trim_end_matches(|c: char| c == 'X' || c == 'x').to_string();
    }

    let n_restrictions =
        (front_restriction.is_some() as u8) + (back_restriction.is_some() as u8);
    if n_restrictions > 1 {
        return Err(ParseError::InvalidSpec(
            "You cannot use multiple placement restrictions for an adapter at the same time. \
             Choose one of ^ADAPTER, ADAPTER$, XADAPTER or ADAPTERX"
                .to_string(),
        ));
    }

    Ok((front_restriction, back_restriction, s))
}

// ---------------------------------------------------------------------------
// Normalize ellipsis (for partial linked-adapter specs)
// ---------------------------------------------------------------------------

fn normalize_ellipsis(
    spec1: &str,
    spec2: &str,
    adapter_type: AdapterType,
) -> Result<(String, AdapterType), ParseError> {
    if adapter_type == AdapterType::Anywhere {
        return Err(ParseError::InvalidSpec(
            "No ellipsis (\"...\") allowed in \"anywhere\" adapters".to_string(),
        ));
    }
    if spec1.is_empty() {
        if adapter_type == AdapterType::Back {
            // -a ...ADAPTER
            Ok((spec2.to_string(), AdapterType::Back))
        } else {
            // -g ...ADAPTER — invalid
            Err(ParseError::InvalidSpec(
                "Invalid adapter specification".to_string(),
            ))
        }
    } else if spec2.is_empty() {
        if adapter_type == AdapterType::Back {
            // -a ADAPTER...
            Ok((spec1.to_string(), AdapterType::Front))
        } else {
            // -g ADAPTER...
            Ok((spec1.to_string(), AdapterType::Front))
        }
    } else {
        Err(ParseError::InvalidSpec(
            "Expected either spec1 or spec2 to be empty".to_string(),
        ))
    }
}

// ---------------------------------------------------------------------------
// Building adapters
// ---------------------------------------------------------------------------

/// Default search parameters used when creating adapters.
#[derive(Debug, Clone)]
pub struct SearchParameters {
    pub max_errors: f64,
    pub min_overlap: usize,
    pub read_wildcards: bool,
    pub adapter_wildcards: bool,
    pub indels: bool,
}

impl Default for SearchParameters {
    fn default() -> Self {
        Self {
            max_errors: 0.1,
            min_overlap: 3,
            read_wildcards: false,
            adapter_wildcards: true,
            indels: true,
        }
    }
}

/// Build `SingleAdapterParams` from an `AdapterSpecification` merged with default search params.
fn build_single_adapter_params(
    sequence: &str,
    name: Option<String>,
    search_params: &SearchParameters,
    overrides: &HashMap<String, ParamValue>,
) -> SingleAdapterParams {
    let mut p = SingleAdapterParams {
        sequence: sequence.to_string(),
        max_errors: search_params.max_errors,
        min_overlap: search_params.min_overlap,
        read_wildcards: search_params.read_wildcards,
        adapter_wildcards: search_params.adapter_wildcards,
        name,
        indels: search_params.indels,
    };

    for (key, value) in overrides {
        match key.as_str() {
            "max_errors" => {
                if let Some(f) = value.as_float() {
                    p.max_errors = f;
                }
            }
            "min_overlap" => {
                if let Some(i) = value.as_int() {
                    p.min_overlap = i as usize;
                }
            }
            "indels" => {
                if let Some(b) = value.as_bool() {
                    p.indels = b;
                }
            }
            _ => {}
        }
    }

    p
}

/// Create a `SingleAdapter` from an `AdapterClass` and params, with optional force_anywhere.
fn create_single_adapter(
    class: AdapterClass,
    params: SingleAdapterParams,
    force_anywhere: bool,
) -> SingleAdapter {
    match class {
        AdapterClass::Front => {
            SingleAdapter::Front(FrontAdapter::with_force_anywhere(params, force_anywhere))
        }
        AdapterClass::Back => {
            SingleAdapter::Back(BackAdapter::with_force_anywhere(params, force_anywhere))
        }
        AdapterClass::Anywhere => SingleAdapter::Anywhere(AnywhereAdapter::new(params)),
        AdapterClass::Prefix => {
            SingleAdapter::Prefix(PrefixAdapter::with_force_anywhere(params, force_anywhere))
        }
        AdapterClass::Suffix => {
            SingleAdapter::Suffix(SuffixAdapter::with_force_anywhere(params, force_anywhere))
        }
        AdapterClass::NonInternalFront => SingleAdapter::NonInternalFront(
            NonInternalFrontAdapter::with_force_anywhere(params, force_anywhere),
        ),
        AdapterClass::NonInternalBack => SingleAdapter::NonInternalBack(
            NonInternalBackAdapter::with_force_anywhere(params, force_anywhere),
        ),
        AdapterClass::RightmostFront => SingleAdapter::RightmostFront(
            RightmostFrontAdapter::with_force_anywhere(params, force_anywhere),
        ),
        AdapterClass::RightmostBack => SingleAdapter::RightmostBack(
            RightmostBackAdapter::with_force_anywhere(params, force_anywhere),
        ),
    }
}

// ---------------------------------------------------------------------------
// make_adapter
// ---------------------------------------------------------------------------

/// Parse an adapter specification (not using `file:` notation) and return
/// an `Adapter`.
pub fn make_adapter(
    spec: &str,
    adapter_type: AdapterType,
    search_params: &SearchParameters,
    name: Option<String>,
) -> Result<Adapter, ParseError> {
    // Check for linked adapter (FRONT...BACK)
    if let Some(dot_pos) = spec.find("...") {
        let spec1 = &spec[..dot_pos];
        let spec2 = &spec[dot_pos + 3..];
        if !spec1.is_empty() && !spec2.is_empty() {
            return make_linked_adapter(spec1, spec2, name, adapter_type, search_params);
        }

        // Partial linked notation: just one side
        let (new_spec, new_type) = normalize_ellipsis(spec1, spec2, adapter_type)?;
        return make_not_linked_adapter(&new_spec, name, new_type, search_params);
    }

    make_not_linked_adapter(spec, name, adapter_type, search_params)
}

/// Create a linked adapter from two specification strings.
fn make_linked_adapter(
    spec1: &str,
    spec2: &str,
    name: Option<String>,
    adapter_type: AdapterType,
    search_params: &SearchParameters,
) -> Result<Adapter, ParseError> {
    if adapter_type == AdapterType::Anywhere {
        return Err(ParseError::InvalidSpec(
            "'anywhere' (-b) adapters may not be linked".to_string(),
        ));
    }

    let front_spec = AdapterSpecification::parse(spec1, AdapterType::Front)?;
    let back_spec = AdapterSpecification::parse(spec2, AdapterType::Back)?;

    let resolved_name = name.or_else(|| front_spec.name.clone());

    let front_anchored = front_spec.restriction.is_some();
    let back_anchored = back_spec.restriction.is_some();

    let mut front_params_map = front_spec.parameters.clone();
    let mut back_params_map = back_spec.parameters.clone();

    let (mut front_required, mut back_required) = if adapter_type == AdapterType::Front {
        (true, true)
    } else {
        (front_anchored, back_anchored)
    };

    // Handle 'required' parameter overrides
    if let Some(v) = front_params_map.remove("required") {
        front_required = v.as_bool().unwrap_or(front_required);
    }
    if let Some(v) = back_params_map.remove("required") {
        back_required = v.as_bool().unwrap_or(back_required);
    }

    let front_adapter_params = build_single_adapter_params(
        &front_spec.sequence,
        Some("linked_front".to_string()),
        search_params,
        &front_params_map,
    );
    let back_adapter_params = build_single_adapter_params(
        &back_spec.sequence,
        Some("linked_back".to_string()),
        search_params,
        &back_params_map,
    );

    let front_adapter = create_single_adapter(front_spec.adapter_class(), front_adapter_params, false);
    let back_adapter = create_single_adapter(back_spec.adapter_class(), back_adapter_params, false);

    let linked = LinkedAdapter::new(
        front_adapter,
        back_adapter,
        front_required,
        back_required,
        resolved_name,
    );

    Ok(Adapter::Linked(linked))
}

/// Create a non-linked adapter from a specification string.
fn make_not_linked_adapter(
    spec: &str,
    name: Option<String>,
    adapter_type: AdapterType,
    search_params: &SearchParameters,
) -> Result<Adapter, ParseError> {
    let aspec = AdapterSpecification::parse(spec, adapter_type)?;
    let adapter_class = aspec.adapter_class();

    let mut params_map = aspec.parameters.clone();

    // Check for `anywhere` parameter overriding to force_anywhere
    let force_anywhere = params_map.remove("anywhere").is_some()
        && matches!(
            adapter_class,
            AdapterClass::Front
                | AdapterClass::Back
                | AdapterClass::RightmostFront
                | AdapterClass::RightmostBack
        );

    if params_map.contains_key("required") {
        return Err(ParseError::InvalidSpec(
            "'required' and 'optional' can only be used within linked adapters".to_string(),
        ));
    }

    let resolved_name = name.or(aspec.name);

    let adapter_params =
        build_single_adapter_params(&aspec.sequence, resolved_name, search_params, &params_map);

    let single = create_single_adapter(adapter_class, adapter_params, force_anywhere);
    Ok(Adapter::Single(single))
}

// ---------------------------------------------------------------------------
// read_adapters_fasta
// ---------------------------------------------------------------------------

/// Read adapter sequences from a FASTA file.
///
/// Returns a `Vec` of `(name, sequence)` pairs.
pub fn read_adapters_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<(Option<String>, String)>, ParseError> {
    let file = std::fs::File::open(path.as_ref())?;
    let reader = std::io::BufReader::new(file);

    let mut records = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            // Save previous record
            if current_name.is_some() || !current_seq.is_empty() {
                records.push((current_name.take(), std::mem::take(&mut current_seq)));
            }
            let header = line[1..].trim().to_string();
            let name = header.split_whitespace().next().map(|s| s.to_string());
            current_name = name;
            current_seq.clear();
        } else {
            current_seq.push_str(&line);
        }
    }
    // Save last record
    if current_name.is_some() || !current_seq.is_empty() {
        records.push((current_name, current_seq));
    }

    Ok(records)
}

// ---------------------------------------------------------------------------
// make_adapters_from_specifications
// ---------------------------------------------------------------------------

/// Create a list of `Adapter` objects from `(adapter_type_str, spec_string)` pairs.
///
/// `type_spec_pairs` is a list of `("front"|"back"|"anywhere", "SPEC")` pairs.
/// `search_params` provides default search parameters that can be overridden by
/// per-adapter parameters in the specification string.
pub fn make_adapters_from_specifications(
    type_spec_pairs: &[(AdapterType, &str)],
    search_params: &SearchParameters,
) -> Result<Vec<Adapter>, ParseError> {
    let mut adapters = Vec::new();
    for &(adapter_type, spec) in type_spec_pairs {
        let mut new_adapters =
            make_adapters_from_one_specification(spec, adapter_type, search_params)?;
        adapters.append(&mut new_adapters);
    }
    Ok(adapters)
}

/// Parse a single adapter specification and return adapter(s).
///
/// Handles `file:path.fasta` notation (including `^file:` and `file$:` prefixes).
fn make_adapters_from_one_specification(
    spec: &str,
    adapter_type: AdapterType,
    search_params: &SearchParameters,
) -> Result<Vec<Adapter>, ParseError> {
    if spec.starts_with("file:")
        || spec.starts_with("^file:")
        || spec.starts_with("file$:")
    {
        let anchoring_prefix;
        let anchoring_suffix;
        let file_spec;

        if spec.starts_with('^') {
            anchoring_prefix = "^";
            anchoring_suffix = "";
            file_spec = &spec[1..]; // "file:..."
        } else if spec.starts_with("file$:") {
            anchoring_prefix = "";
            anchoring_suffix = "$";
            let (path_part, _, params_spec) = spec[6..].partition_at_char(';');
            let records = read_adapters_fasta(path_part)?;
            let search_overrides = parse_search_parameters(params_spec)?;
            let mut merged_params = search_params.clone();
            apply_overrides(&mut merged_params, &search_overrides);
            let mut adapters = Vec::new();
            for (name, seq) in records {
                let full_spec = format!("{}{}{}", anchoring_prefix, seq, anchoring_suffix);
                let adapter = make_adapter(&full_spec, adapter_type, &merged_params, name)?;
                adapters.push(adapter);
            }
            return Ok(adapters);
        } else {
            anchoring_prefix = "";
            anchoring_suffix = "";
            file_spec = spec;
        }

        // file_spec starts with "file:"
        let after_file = &file_spec[5..];
        let (path_part, _, params_spec) = after_file.partition_at_char(';');
        let records = read_adapters_fasta(path_part)?;
        let search_overrides = parse_search_parameters(params_spec)?;
        let mut merged_params = search_params.clone();
        apply_overrides(&mut merged_params, &search_overrides);
        let mut adapters = Vec::new();
        for (name, seq) in records {
            let full_spec = format!("{}{}{}", anchoring_prefix, seq, anchoring_suffix);
            let adapter = make_adapter(&full_spec, adapter_type, &merged_params, name)?;
            adapters.push(adapter);
        }
        Ok(adapters)
    } else {
        let adapter = make_adapter(spec, adapter_type, search_params, None)?;
        Ok(vec![adapter])
    }
}

fn apply_overrides(params: &mut SearchParameters, overrides: &HashMap<String, ParamValue>) {
    for (key, value) in overrides {
        match key.as_str() {
            "max_errors" => {
                if let Some(f) = value.as_float() {
                    params.max_errors = f;
                }
            }
            "min_overlap" => {
                if let Some(i) = value.as_int() {
                    params.min_overlap = i as usize;
                }
            }
            "indels" => {
                if let Some(b) = value.as_bool() {
                    params.indels = b;
                }
            }
            _ => {}
        }
    }
}

/// Helper trait for partitioning a string at the first occurrence of a character.
trait PartitionAtChar {
    fn partition_at_char(&self, ch: char) -> (&str, &str, &str);
}

impl PartitionAtChar for str {
    fn partition_at_char(&self, ch: char) -> (&str, &str, &str) {
        if let Some(pos) = self.find(ch) {
            (&self[..pos], &self[pos..pos + 1], &self[pos + 1..])
        } else {
            (self, "", "")
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- expand_braces tests --

    #[test]
    fn test_expand_braces_basic() {
        assert_eq!(expand_braces("TGA{5}CT").unwrap(), "TGAAAAACT");
    }

    #[test]
    fn test_expand_braces_no_braces() {
        assert_eq!(expand_braces("ACGT").unwrap(), "ACGT");
    }

    #[test]
    fn test_expand_braces_zero() {
        assert_eq!(expand_braces("TGA{0}CT").unwrap(), "TGCT");
    }

    #[test]
    fn test_expand_braces_one() {
        assert_eq!(expand_braces("TGA{1}CT").unwrap(), "TGACT");
    }

    #[test]
    fn test_expand_braces_multiple() {
        assert_eq!(expand_braces("A{3}T{2}").unwrap(), "AAATT");
    }

    #[test]
    fn test_expand_braces_error_leading_brace() {
        assert!(expand_braces("{5}CT").is_err());
    }

    #[test]
    fn test_expand_braces_error_unmatched_close() {
        assert!(expand_braces("A}CT").is_err());
    }

    #[test]
    fn test_expand_braces_error_unterminated() {
        assert!(expand_braces("A{5").is_err());
    }

    #[test]
    fn test_expand_braces_error_too_large() {
        assert!(expand_braces("A{10001}").is_err());
    }

    #[test]
    fn test_expand_braces_empty() {
        assert_eq!(expand_braces("").unwrap(), "");
    }

    // -- parse_search_parameters tests --

    #[test]
    fn test_parse_search_params_empty() {
        let result = parse_search_parameters("").unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_parse_search_params_basic() {
        let result = parse_search_parameters("e=0.1;o=3").unwrap();
        assert_eq!(result.get("max_errors"), Some(&ParamValue::Float(0.1)));
        assert_eq!(result.get("min_overlap"), Some(&ParamValue::Int(3)));
    }

    #[test]
    fn test_parse_search_params_bool() {
        let result = parse_search_parameters("anywhere").unwrap();
        assert_eq!(result.get("anywhere"), Some(&ParamValue::Bool(true)));
    }

    #[test]
    fn test_parse_search_params_noindels_normalization() {
        let result = parse_search_parameters("noindels").unwrap();
        assert_eq!(result.get("indels"), Some(&ParamValue::Bool(false)));
        assert!(!result.contains_key("noindels"));
    }

    #[test]
    fn test_parse_search_params_optional_normalization() {
        let result = parse_search_parameters("optional").unwrap();
        assert_eq!(result.get("required"), Some(&ParamValue::Bool(false)));
        assert!(!result.contains_key("optional"));
    }

    #[test]
    fn test_parse_search_params_unknown_key() {
        assert!(parse_search_parameters("unknown=1").is_err());
    }

    #[test]
    fn test_parse_search_params_duplicate_key() {
        assert!(parse_search_parameters("e=0.1;e=0.2").is_err());
    }

    #[test]
    fn test_parse_search_params_optional_and_required() {
        assert!(parse_search_parameters("optional;required").is_err());
    }

    #[test]
    fn test_parse_search_params_indels_and_noindels() {
        assert!(parse_search_parameters("indels;noindels").is_err());
    }

    #[test]
    fn test_parse_search_params_missing_value() {
        assert!(parse_search_parameters("e=").is_err());
    }

    // -- AdapterSpecification::parse tests --

    #[test]
    fn test_parse_basic_back() {
        let spec = AdapterSpecification::parse("ACGT", AdapterType::Back).unwrap();
        assert_eq!(spec.name, None);
        assert_eq!(spec.restriction, None);
        assert_eq!(spec.sequence, "ACGT");
        assert_eq!(spec.adapter_type, AdapterType::Back);
        assert!(!spec.rightmost);
    }

    #[test]
    fn test_parse_named() {
        let spec = AdapterSpecification::parse("myname=ACGT", AdapterType::Back).unwrap();
        assert_eq!(spec.name, Some("myname".to_string()));
        assert_eq!(spec.sequence, "ACGT");
    }

    #[test]
    fn test_parse_with_params() {
        let spec =
            AdapterSpecification::parse("a_name=ACGT;anywhere", AdapterType::Back).unwrap();
        assert_eq!(spec.name, Some("a_name".to_string()));
        assert_eq!(spec.sequence, "ACGT");
        assert_eq!(
            spec.parameters.get("anywhere"),
            Some(&ParamValue::Bool(true))
        );
    }

    #[test]
    fn test_parse_anchored_front() {
        let spec = AdapterSpecification::parse("^ACGT", AdapterType::Front).unwrap();
        assert_eq!(spec.restriction, Some(Restriction::Anchored));
        assert_eq!(spec.sequence, "ACGT");
    }

    #[test]
    fn test_parse_anchored_back() {
        let spec = AdapterSpecification::parse("ACGT$", AdapterType::Back).unwrap();
        assert_eq!(spec.restriction, Some(Restriction::Anchored));
        assert_eq!(spec.sequence, "ACGT");
    }

    #[test]
    fn test_parse_noninternal_front() {
        let spec = AdapterSpecification::parse("XACGT", AdapterType::Front).unwrap();
        assert_eq!(spec.restriction, Some(Restriction::NonInternal));
        assert_eq!(spec.sequence, "ACGT");
    }

    #[test]
    fn test_parse_noninternal_back() {
        let spec = AdapterSpecification::parse("ACGTX", AdapterType::Back).unwrap();
        assert_eq!(spec.restriction, Some(Restriction::NonInternal));
        assert_eq!(spec.sequence, "ACGT");
    }

    #[test]
    fn test_parse_front_back_restriction_error() {
        // Can't use $ with front adapter
        assert!(AdapterSpecification::parse("ACGT$", AdapterType::Front).is_err());
    }

    #[test]
    fn test_parse_back_front_restriction_error() {
        // Can't use ^ with back adapter
        assert!(AdapterSpecification::parse("^ACGT", AdapterType::Back).is_err());
    }

    #[test]
    fn test_parse_anywhere_restriction_error() {
        assert!(AdapterSpecification::parse("^ACGT", AdapterType::Anywhere).is_err());
    }

    #[test]
    fn test_parse_only_x() {
        // Special case: only X characters
        let spec = AdapterSpecification::parse("XXX", AdapterType::Back).unwrap();
        assert_eq!(spec.sequence, "XXX");
        assert_eq!(spec.restriction, None);
    }

    #[test]
    fn test_parse_brace_expansion() {
        let spec = AdapterSpecification::parse("A{3}T", AdapterType::Back).unwrap();
        assert_eq!(spec.sequence, "AAAT");
    }

    #[test]
    fn test_parse_rightmost() {
        let spec =
            AdapterSpecification::parse("ACGT;rightmost", AdapterType::Front).unwrap();
        assert!(spec.rightmost);
        assert_eq!(spec.adapter_class(), AdapterClass::RightmostFront);
    }

    #[test]
    fn test_parse_rightmost_with_restriction_error() {
        assert!(
            AdapterSpecification::parse("^ACGT;rightmost", AdapterType::Front).is_err()
        );
    }

    #[test]
    fn test_parse_min_overlap_anchored_error() {
        assert!(
            AdapterSpecification::parse("^ACGT;o=2", AdapterType::Front).is_err()
        );
    }

    #[test]
    fn test_parse_min_overlap_clamped() {
        let spec = AdapterSpecification::parse("ACG;o=100", AdapterType::Back).unwrap();
        // min_overlap should be clamped to sequence length (3)
        assert_eq!(
            spec.parameters.get("min_overlap"),
            Some(&ParamValue::Int(3))
        );
    }

    // -- adapter_class tests --

    #[test]
    fn test_adapter_class_front() {
        let spec = AdapterSpecification::parse("ACGT", AdapterType::Front).unwrap();
        assert_eq!(spec.adapter_class(), AdapterClass::Front);
    }

    #[test]
    fn test_adapter_class_back() {
        let spec = AdapterSpecification::parse("ACGT", AdapterType::Back).unwrap();
        assert_eq!(spec.adapter_class(), AdapterClass::Back);
    }

    #[test]
    fn test_adapter_class_anywhere() {
        let spec = AdapterSpecification::parse("ACGT", AdapterType::Anywhere).unwrap();
        assert_eq!(spec.adapter_class(), AdapterClass::Anywhere);
    }

    #[test]
    fn test_adapter_class_prefix() {
        let spec = AdapterSpecification::parse("^ACGT", AdapterType::Front).unwrap();
        assert_eq!(spec.adapter_class(), AdapterClass::Prefix);
    }

    #[test]
    fn test_adapter_class_suffix() {
        let spec = AdapterSpecification::parse("ACGT$", AdapterType::Back).unwrap();
        assert_eq!(spec.adapter_class(), AdapterClass::Suffix);
    }

    #[test]
    fn test_adapter_class_noninternal_front() {
        let spec = AdapterSpecification::parse("XACGT", AdapterType::Front).unwrap();
        assert_eq!(spec.adapter_class(), AdapterClass::NonInternalFront);
    }

    #[test]
    fn test_adapter_class_noninternal_back() {
        let spec = AdapterSpecification::parse("ACGTX", AdapterType::Back).unwrap();
        assert_eq!(spec.adapter_class(), AdapterClass::NonInternalBack);
    }

    // -- make_adapter tests --

    #[test]
    fn test_make_adapter_simple_back() {
        let params = SearchParameters::default();
        let adapter = make_adapter("AGATCGGAAGAGC", AdapterType::Back, &params, None).unwrap();
        match &adapter {
            Adapter::Single(SingleAdapter::Back(_)) => {}
            _ => panic!("Expected BackAdapter"),
        }
    }

    #[test]
    fn test_make_adapter_simple_front() {
        let params = SearchParameters::default();
        let adapter = make_adapter("AGATCGGAAGAGC", AdapterType::Front, &params, None).unwrap();
        match &adapter {
            Adapter::Single(SingleAdapter::Front(_)) => {}
            _ => panic!("Expected FrontAdapter"),
        }
    }

    #[test]
    fn test_make_adapter_anchored_front() {
        let params = SearchParameters::default();
        let adapter = make_adapter("^ACGT", AdapterType::Front, &params, None).unwrap();
        match &adapter {
            Adapter::Single(SingleAdapter::Prefix(_)) => {}
            _ => panic!("Expected PrefixAdapter"),
        }
    }

    #[test]
    fn test_make_adapter_anchored_back() {
        let params = SearchParameters::default();
        let adapter = make_adapter("ACGT$", AdapterType::Back, &params, None).unwrap();
        match &adapter {
            Adapter::Single(SingleAdapter::Suffix(_)) => {}
            _ => panic!("Expected SuffixAdapter"),
        }
    }

    #[test]
    fn test_make_adapter_noninternal_front() {
        let params = SearchParameters::default();
        let adapter = make_adapter("XACGT", AdapterType::Front, &params, None).unwrap();
        match &adapter {
            Adapter::Single(SingleAdapter::NonInternalFront(_)) => {}
            _ => panic!("Expected NonInternalFrontAdapter"),
        }
    }

    #[test]
    fn test_make_adapter_noninternal_back() {
        let params = SearchParameters::default();
        let adapter = make_adapter("ACGTX", AdapterType::Back, &params, None).unwrap();
        match &adapter {
            Adapter::Single(SingleAdapter::NonInternalBack(_)) => {}
            _ => panic!("Expected NonInternalBackAdapter"),
        }
    }

    #[test]
    fn test_make_adapter_anywhere() {
        let params = SearchParameters::default();
        let adapter = make_adapter("ACGT", AdapterType::Anywhere, &params, None).unwrap();
        match &adapter {
            Adapter::Single(SingleAdapter::Anywhere(_)) => {}
            _ => panic!("Expected AnywhereAdapter"),
        }
    }

    #[test]
    fn test_make_adapter_rightmost_front() {
        let params = SearchParameters::default();
        let adapter =
            make_adapter("ACGT;rightmost", AdapterType::Front, &params, None).unwrap();
        match &adapter {
            Adapter::Single(SingleAdapter::RightmostFront(_)) => {}
            _ => panic!("Expected RightmostFrontAdapter"),
        }
    }

    #[test]
    fn test_make_adapter_rightmost_back() {
        let params = SearchParameters::default();
        let adapter =
            make_adapter("ACGT;rightmost", AdapterType::Back, &params, None).unwrap();
        match &adapter {
            Adapter::Single(SingleAdapter::RightmostBack(_)) => {}
            _ => panic!("Expected RightmostBackAdapter"),
        }
    }

    #[test]
    fn test_make_adapter_linked() {
        let params = SearchParameters::default();
        let adapter =
            make_adapter("AAAA...TTTT", AdapterType::Back, &params, None).unwrap();
        match &adapter {
            Adapter::Linked(_) => {}
            _ => panic!("Expected LinkedAdapter"),
        }
    }

    #[test]
    fn test_make_adapter_linked_front_type() {
        let params = SearchParameters::default();
        let adapter =
            make_adapter("AAAA...TTTT", AdapterType::Front, &params, None).unwrap();
        match &adapter {
            Adapter::Linked(_) => {}
            _ => panic!("Expected LinkedAdapter"),
        }
    }

    #[test]
    fn test_make_adapter_linked_anywhere_error() {
        let params = SearchParameters::default();
        let result = make_adapter("AAAA...TTTT", AdapterType::Anywhere, &params, None);
        assert!(result.is_err());
    }

    #[test]
    fn test_make_adapter_ellipsis_back() {
        // -a ...ADAPTER -> BackAdapter
        let params = SearchParameters::default();
        let adapter = make_adapter("...TTTT", AdapterType::Back, &params, None).unwrap();
        match &adapter {
            Adapter::Single(SingleAdapter::Back(_)) => {}
            _ => panic!("Expected BackAdapter for ...ADAPTER with -a"),
        }
    }

    #[test]
    fn test_make_adapter_ellipsis_front() {
        // -a ADAPTER... -> FrontAdapter
        let params = SearchParameters::default();
        let adapter = make_adapter("AAAA...", AdapterType::Back, &params, None).unwrap();
        match &adapter {
            Adapter::Single(SingleAdapter::Front(_)) => {}
            _ => panic!("Expected FrontAdapter for ADAPTER... with -a"),
        }
    }

    #[test]
    fn test_make_adapter_with_name() {
        let params = SearchParameters::default();
        let adapter =
            make_adapter("my_adapter=ACGT", AdapterType::Back, &params, None).unwrap();
        assert_eq!(adapter.name(), "my_adapter");
    }

    #[test]
    fn test_make_adapter_external_name() {
        let params = SearchParameters::default();
        let adapter = make_adapter(
            "ACGT",
            AdapterType::Back,
            &params,
            Some("ext_name".to_string()),
        )
        .unwrap();
        assert_eq!(adapter.name(), "ext_name");
    }

    #[test]
    fn test_make_adapter_force_anywhere() {
        let params = SearchParameters::default();
        let adapter =
            make_adapter("ACGT;anywhere", AdapterType::Back, &params, None).unwrap();
        // Should still be a BackAdapter, but with force_anywhere enabled internally
        match &adapter {
            Adapter::Single(SingleAdapter::Back(_)) => {}
            _ => panic!("Expected BackAdapter with force_anywhere"),
        }
    }

    // -- make_adapters_from_specifications tests --

    #[test]
    fn test_make_adapters_from_specs() {
        let params = SearchParameters::default();
        let specs = vec![
            (AdapterType::Back, "AGATCGGAAGAGC"),
            (AdapterType::Front, "TGGAATTCTCGG"),
        ];
        let adapters = make_adapters_from_specifications(&specs, &params).unwrap();
        assert_eq!(adapters.len(), 2);
    }

    // -- read_adapters_fasta tests --

    #[test]
    fn test_read_adapters_fasta() {
        use std::io::Write;
        let dir = std::env::temp_dir().join("cutadapt_test_fasta");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test_adapters.fasta");
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, ">adapter1 some description").unwrap();
            writeln!(f, "AGATCGGAAGAGC").unwrap();
            writeln!(f, ">adapter2").unwrap();
            writeln!(f, "TGGAATTCTCGG").unwrap();
        }

        let records = read_adapters_fasta(&path).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].0, Some("adapter1".to_string()));
        assert_eq!(records[0].1, "AGATCGGAAGAGC");
        assert_eq!(records[1].0, Some("adapter2".to_string()));
        assert_eq!(records[1].1, "TGGAATTCTCGG");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_read_adapters_fasta_multiline() {
        use std::io::Write;
        let dir = std::env::temp_dir().join("cutadapt_test_fasta_ml");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test_multiline.fasta");
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, ">adapter1").unwrap();
            writeln!(f, "AGATC").unwrap();
            writeln!(f, "GGAAG").unwrap();
            writeln!(f, ">adapter2").unwrap();
            writeln!(f, "TGGAA").unwrap();
        }

        let records = read_adapters_fasta(&path).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].1, "AGATCGGAAG");
        assert_eq!(records[1].1, "TGGAA");

        std::fs::remove_dir_all(&dir).ok();
    }

    // -- restriction parsing tests --

    #[test]
    fn test_parse_restrictions_none() {
        let (f, b, s) = parse_restrictions("ACGT").unwrap();
        assert_eq!(f, None);
        assert_eq!(b, None);
        assert_eq!(s, "ACGT");
    }

    #[test]
    fn test_parse_restrictions_front_anchored() {
        let (f, b, s) = parse_restrictions("^ACGT").unwrap();
        assert_eq!(f, Some(Restriction::Anchored));
        assert_eq!(b, None);
        assert_eq!(s, "ACGT");
    }

    #[test]
    fn test_parse_restrictions_back_anchored() {
        let (f, b, s) = parse_restrictions("ACGT$").unwrap();
        assert_eq!(f, None);
        assert_eq!(b, Some(Restriction::Anchored));
        assert_eq!(s, "ACGT");
    }

    #[test]
    fn test_parse_restrictions_front_noninternal() {
        let (f, b, s) = parse_restrictions("XACGT").unwrap();
        assert_eq!(f, Some(Restriction::NonInternal));
        assert_eq!(b, None);
        assert_eq!(s, "ACGT");
    }

    #[test]
    fn test_parse_restrictions_back_noninternal() {
        let (f, b, s) = parse_restrictions("ACGTX").unwrap();
        assert_eq!(f, None);
        assert_eq!(b, Some(Restriction::NonInternal));
        assert_eq!(s, "ACGT");
    }

    #[test]
    fn test_parse_restrictions_both_error() {
        assert!(parse_restrictions("^ACGT$").is_err());
    }

    #[test]
    fn test_parse_restrictions_front_double_error() {
        // ^X would be two front restrictions
        assert!(parse_restrictions("^XACGT").is_err());
    }

    // -- normalize_ellipsis tests --

    #[test]
    fn test_normalize_ellipsis_back_left_empty() {
        let (s, t) = normalize_ellipsis("", "ADAPTER", AdapterType::Back).unwrap();
        assert_eq!(s, "ADAPTER");
        assert_eq!(t, AdapterType::Back);
    }

    #[test]
    fn test_normalize_ellipsis_back_right_empty() {
        let (s, t) = normalize_ellipsis("ADAPTER", "", AdapterType::Back).unwrap();
        assert_eq!(s, "ADAPTER");
        assert_eq!(t, AdapterType::Front);
    }

    #[test]
    fn test_normalize_ellipsis_front_right_empty() {
        let (s, t) = normalize_ellipsis("ADAPTER", "", AdapterType::Front).unwrap();
        assert_eq!(s, "ADAPTER");
        assert_eq!(t, AdapterType::Front);
    }

    #[test]
    fn test_normalize_ellipsis_front_left_empty_error() {
        assert!(normalize_ellipsis("", "ADAPTER", AdapterType::Front).is_err());
    }

    #[test]
    fn test_normalize_ellipsis_anywhere_error() {
        assert!(normalize_ellipsis("A", "B", AdapterType::Anywhere).is_err());
    }
}
